import subprocess
import sys
from pathlib import Path
from collections import namedtuple

from jug import TaskGenerator
from jug.utils import timed_path


WORK_DIR = Path('.').resolve()
PROJECT_ROOT = WORK_DIR.parent.parent
THREADS = "8"
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

for sample in SAMPLES:
    paths = sample_paths(sample)
    long_bam = align_long_reads(sample, REFERENCE, paths.ont)
    short_bam = align_short_reads(sample, REFERENCE, paths.ilm1, paths.ilm2)
    index_bam(sort_bam(long_bam))
    index_bam(sort_bam(short_bam))
