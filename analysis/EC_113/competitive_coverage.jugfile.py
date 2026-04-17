import subprocess
from pathlib import Path

from jug import TaskGenerator
from jug.utils import timed_path

from ncpus import get_ncpus

WORK_DIR = Path('.').resolve()
PROJECT_ROOT = WORK_DIR.parent.parent
THREADS = str(get_ncpus())

MAPPED_SP_DIR = PROJECT_ROOT / "resource_generation" / "ShortRead_mappings" / "outputs" / "mapped_sp"
FILTER_NGL = str(WORK_DIR / "filter_SHD1_0457_1.ngl")


@TaskGenerator
def filter_bam(sample, input_bam, ngl_script):
    """Filter competitively-mapped BAM to keep only reads mapped to SHD1_0457_1."""
    output_dir = WORK_DIR / "outputs" / "mapped_sp_filtered"
    output_dir.mkdir(exist_ok=True, parents=True)

    output = output_dir / f"{sample}_ShanghaiDogsMAGsSpecies.filtered_SHD1_0457_1.sorted.bam"
    subprocess.run(
        ["ngless", f"-j{THREADS}", ngl_script, str(input_bam), str(output)],
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
    output_dir = WORK_DIR / "outputs" / "coverage_competitive"
    output_dir.mkdir(exist_ok=True, parents=True)
    output = output_dir / f"{sorted_bam.stem}.bedgraph"
    with open(output, "w") as f:
        subprocess.run(
            ["bedtools", "genomecov", "-ibam", str(sorted_bam), "-bga"],
            stdout=f,
            check=True,
        )
    return str(output)


@TaskGenerator
def plot_coverage_profiles(bedgraphs, samples):
    """Plot coverage profiles for competitively-mapped short reads across the reference."""
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt

    WINDOW_SIZE = 1000

    def bedgraph_to_binned(path):
        regions = []
        max_end = 0
        with open(path) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if parts[0] != "SHD1_0457_1":
                    continue
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

    fig, ax = plt.subplots(figsize=(14, 5))
    for sample, bg_path in zip(samples, bedgraphs):
        positions, coverage = bedgraph_to_binned(bg_path)
        ax.plot(positions, coverage, alpha=0.5, linewidth=0.5, label=sample)
    ax.set_ylabel("Mean coverage")
    ax.set_xlabel("Position (Mb)")
    ax.set_title("Coverage profiles across SHD1_0457 (competitive mapping)")
    ax.legend(fontsize=7, ncol=3, loc="upper right")
    fig.tight_layout()

    output_path = output_dir / "coverage_profiles_competitive.png"
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return str(output_path)


SAMPLES = [
     'D003',
     'D004',
     'D005',
     'D006',
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

bedgraphs = []
for sample in SAMPLES:
    input_bam = timed_path(str(MAPPED_SP_DIR / f"{sample}_ShanghaiDogsMAGsSpecies.bam"))
    filtered = filter_bam(sample, input_bam, FILTER_NGL)
    index_bam(filtered)
    bedgraphs.append(bedtools_genomecov(filtered))

plot_coverage_profiles(bedgraphs, SAMPLES)
