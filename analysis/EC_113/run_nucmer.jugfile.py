import re
import subprocess
from itertools import combinations
from pathlib import Path

from jug import TaskGenerator
from jug.utils import timed_path


DATA_DIR = Path("data")
OUTPUTS_DIR = Path("outputs")
LOGS_DIR = Path("logs")
OUTPUT_DIR = OUTPUTS_DIR / "nucmer"
ID_RE = re.compile(r"^SHD1_(\d+)\.fna$")
ID_ECE_RE = re.compile(r"^SHD1_EC\.(\d+)\.fna$")


def extract_id(path):
    if match := ID_RE.match(path.name):
        return match.group(1)
    if match := ID_ECE_RE.match(path.name):
        n = match.group(1)
        return f"ECE.{n}"
    raise ValueError(f"Could not extract ID from filename: {path.name}")

def log_paths_for(output_path):
    relative_path = Path(output_path).relative_to(OUTPUTS_DIR)
    log_base = LOGS_DIR / relative_path
    log_base.parent.mkdir(parents=True, exist_ok=True)
    return (
        log_base.with_name(log_base.name + ".stdout"),
        log_base.with_name(log_base.name + ".stderr"),
    )


def run_logged(command, output_path, **kwargs):
    stdout_path, stderr_path = log_paths_for(output_path)
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        subprocess.run(command, check=True, stdout=stdout, stderr=stderr, **kwargs)


@TaskGenerator
def run_nucmer(file1, file2, id1, id2):
    output_dir = OUTPUTS_DIR / "nucmer"
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = output_dir / f"{id1}_vs_{id2}"
    run_logged(["nucmer", "-p", str(prefix), str(file1), str(file2)], prefix)

    return str(prefix) + ".delta"


@TaskGenerator
def run_dnadiff(delta):
    output_dir = OUTPUTS_DIR / "dnadiff"
    output_dir.mkdir(parents=True, exist_ok=True)

    prefix = Path(delta).stem
    prefix = output_dir / prefix
    run_logged(["dnadiff", "-d", str(delta), "-p", str(prefix)], prefix)
    return str(prefix)


@TaskGenerator
def run_mummerplot(delta):
    output_dir = OUTPUTS_DIR / "mummerplot"
    output_dir.mkdir(parents=True, exist_ok=True)

    prefix = output_dir / Path(delta).stem
    run_logged(["mummerplot", "--png", "-p", str(prefix), str(delta)], prefix)
    return str(prefix) + ".png"


samples = sorted(
    ((extract_id(path), timed_path(str(path))) for path in DATA_DIR.glob("*.fna")),
    key=lambda sample: sample[0],
)

for (id1, file1), (id2, file2) in combinations(samples, 2):
    delta = run_nucmer(file1, file2, id1, id2)
    run_dnadiff(delta)
    run_mummerplot(delta)
