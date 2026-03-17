import subprocess
import sys
from pathlib import Path
from collections import namedtuple

from jug import TaskGenerator
from jug.utils import timed_path

from ncpus import get_ncpus


WORK_DIR = Path('.').resolve()
PROJECT_ROOT = WORK_DIR.parent.parent
THREADS = str(get_ncpus())
REFERENCE = timed_path(str(WORK_DIR / "data" / "SHD1_0457.fna"))


SamplePaths = namedtuple("SamplePaths", ["ont", "ilm1", "ilm2"])

def sample_paths(sample):
    return SamplePaths(
        ont=timed_path(
            str(
                PROJECT_ROOT
                / "data"
                / "ShanghaiDogsFastQ"
                / "ONT"
                / sample
                / f"{sample}_pass.fq.gz"
            )
        ),
        ilm1=timed_path(
            str(
                PROJECT_ROOT
                / "data"
                / "ShanghaiDogsFastQ"
                / "ILM"
                / sample
                / f"{sample}.pair.1.fq.gz"
            )
        ),
        ilm2=timed_path(
            str(
                PROJECT_ROOT
                / "data"
                / "ShanghaiDogsFastQ"
                / "ILM"
                / sample
                / f"{sample}.pair.2.fq.gz"
            )
        ),
    )


def _run_minimap2_to_bam(minimap2_args, output_path):
    """Run minimap2 and pipe output through samtools view to produce BAM."""
    minimap2 = subprocess.Popen(
        minimap2_args,
        stdout=subprocess.PIPE,
    )
    with open(output_path, "wb") as bam:
        samtools = subprocess.Popen(
            ["samtools", "view", "-b", "-@", THREADS],
            stdin=minimap2.stdout,
            stdout=bam,
        )
    minimap2.stdout.close()
    samtools_rc = samtools.wait()
    minimap2_rc = minimap2.wait()
    if minimap2_rc != 0:
        raise subprocess.CalledProcessError(minimap2_rc, minimap2_args)
    if samtools_rc != 0:
        raise subprocess.CalledProcessError(samtools_rc, ["samtools", "view"])


@TaskGenerator
def align_long_reads(sample, reference, ont_reads):
    output_dir = WORK_DIR / "outputs" / "mapped"
    output_dir.mkdir(exist_ok=True)

    output = output_dir / f"{sample}_SHD1_0457_LR.bam"
    _run_minimap2_to_bam(
        [
            "minimap2",
            "-ax",
            "map-ont",
            "-t",
            THREADS,
            str(reference),
            str(ont_reads),
        ],
        output,
    )
    return str(output)


@TaskGenerator
def align_short_reads(sample, reference, ilm_read1, ilm_read2):
    output_dir = WORK_DIR / "outputs" / "mapped"
    output_dir.mkdir(exist_ok=True)

    output = output_dir / f"{sample}_SHD1_0457_SR.bam"
    _run_minimap2_to_bam(
        [
            "minimap2",
            "-ax",
            "sr",
            "-t",
            THREADS,
            str(reference),
            str(ilm_read1),
            str(ilm_read2),
        ],
        output,
    )
    return str(output)


@TaskGenerator
def sort_bam(bam_path):
    bam_path = Path(bam_path)
    output = bam_path.with_suffix(".sorted.bam")
    subprocess.run(
        ["samtools", "sort", "-@", THREADS, "-o", str(output), str(bam_path)],
        check=True,
    )
    return str(output)


@TaskGenerator
def index_bam(bam_path):
    subprocess.run(["samtools", "index", str(bam_path)], check=True)
    return f"{bam_path}.bai"


@TaskGenerator
def bedtools_genomecov(sorted_bam):
    """Run bedtools genomecov to produce a BedGraph coverage file."""
    sorted_bam = Path(sorted_bam)
    output_dir = WORK_DIR / "outputs" / "coverage"
    output_dir.mkdir(exist_ok=True)
    output = output_dir / f"{sorted_bam.stem}.bedgraph"
    with open(output, "w") as f:
        subprocess.run(
            ["bedtools", "genomecov", "-ibam", str(sorted_bam), "-bga"],
            stdout=f,
            check=True,
        )
    return str(output)


@TaskGenerator
def plot_coverage_profiles(lr_bedgraphs, sr_bedgraphs, samples):
    """Plot coverage profiles for long and short reads across the reference."""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    WINDOW_SIZE = 1000

    def bedgraph_to_binned(path):
        regions = []
        max_end = 0
        with open(path) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                start, end, depth = int(parts[1]), int(parts[2]), float(parts[3])
                regions.append((start, end, depth))
                if end > max_end:
                    max_end = end
        cov = np.zeros(max_end, dtype=np.float32)
        for start, end, depth in regions:
            cov[start:end] = depth
        n_bins = max_end // WINDOW_SIZE
        binned = cov[: n_bins * WINDOW_SIZE].reshape(n_bins, WINDOW_SIZE).mean(axis=1)
        positions = np.arange(n_bins) * WINDOW_SIZE / 1e6
        return positions, binned

    output_dir = WORK_DIR / "outputs" / "figures"
    output_dir.mkdir(exist_ok=True, parents=True)

    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    for ax, bedgraphs, label in [
        (axes[0], lr_bedgraphs, "Long reads (ONT)"),
        (axes[1], sr_bedgraphs, "Short reads (Illumina)"),
    ]:
        for sample, bg_path in zip(samples, bedgraphs):
            positions, coverage = bedgraph_to_binned(bg_path)
            ax.plot(positions, coverage, alpha=0.5, linewidth=0.5, label=sample)
        ax.set_ylabel("Mean coverage")
        ax.set_title(label)
        ax.legend(fontsize=7, ncol=3, loc="upper right")

    axes[1].set_xlabel("Position (Mb)")
    fig.suptitle("Coverage profiles across SHD1_0457", fontsize=14)
    fig.tight_layout()

    output_path = output_dir / "coverage_profiles.png"
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return str(output_path)


SAMPLES = [
     'D003',
     'D004',
     'D005',
    # 'D006', # This sample contains two read files, so needs special handling
     'D011',
     'D012',
     'D013',
     'D019',
     'D020',
     'D022',
     'D023',
     'D028',
     'D030',
     'D043',
     'D044',
     ]

lr_bedgraphs = []
sr_bedgraphs = []
for sample in SAMPLES:
    paths = sample_paths(sample)
    long_bam = align_long_reads(sample, REFERENCE, paths.ont)
    short_bam = align_short_reads(sample, REFERENCE, paths.ilm1, paths.ilm2)
    sorted_long = sort_bam(long_bam)
    sorted_short = sort_bam(short_bam)
    index_bam(sorted_long)
    index_bam(sorted_short)
    lr_bedgraphs.append(bedtools_genomecov(sorted_long))
    sr_bedgraphs.append(bedtools_genomecov(sorted_short))

plot_coverage_profiles(lr_bedgraphs, sr_bedgraphs, SAMPLES)
