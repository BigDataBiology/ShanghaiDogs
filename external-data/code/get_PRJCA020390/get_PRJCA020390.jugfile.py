import csv
import subprocess
from pathlib import Path

from jug import TaskGenerator

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR.parent.parent / 'data' / 'PRJCA020390'


@TaskGenerator
def download_file(url, output_path):
    """Download a single file via wget, skipping if it already exists."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        return str(output_path)
    subprocess.run(
        ['wget', '-q', '-O', str(output_path), url],
        check=True,
    )
    return str(output_path)


with open(BASE_DIR / 'RunInfo.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        sample = row['Title']
        urls = row['Download_path'].split('|')
        filenames = row['FileName'].split('|')
        for url, fname in zip(urls, filenames):
            opath = DATA_DIR / sample / fname
            download_file(url, str(opath))
