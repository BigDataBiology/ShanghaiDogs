import subprocess
import sys
from pathlib import Path
from collections import namedtuple

from jug import TaskGenerator
from jug.utils import timed_path


WORK_DIR = Path(__file__).resolve().parent
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


@TaskGenerator
def align_long_reads(sample, reference, ont_reads):
    output = WORK_DIR / "outputs" / "mapped" / f"{sample}_SHD1_0457_LR.sam"
    with output.open("wb") as sam:
        subprocess.run(
            [
                "minimap2",
                "-ax",
                "map-ont",
                "-t",
                THREADS,
                str(reference),
                str(ont_reads),
            ],
            check=True,
            stdout=sam,
        )
    return str(output)


@TaskGenerator
def align_short_reads(sample, reference, ilm_read1, ilm_read2):
    output = WORK_DIR / f"{sample}_SHD1_0457_SR.sam"
    with output.open("wb") as sam:
        subprocess.run(
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
            check=True,
            stdout=sam,
        )
    return str(output)


@TaskGenerator
def sort_bam(sam_path):
    sam_path = Path(sam_path)
    output = sam_path.with_suffix(".sorted.bam")
    with output.open("wb") as bam:
        subprocess.run(
            ["samtools", "sort", "-@", THREADS, str(sam_path)],
            check=True,
            stdout=bam,
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
    long_sam = align_long_reads(sample, REFERENCE, paths.ont)
    short_sam = align_short_reads(sample, REFERENCE, paths.ilm1, paths.ilm2)
    index_bam(sort_bam(long_sam))
    index_bam(sort_bam(short_sam))
