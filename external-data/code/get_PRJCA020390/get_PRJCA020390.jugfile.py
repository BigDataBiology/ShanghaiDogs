import csv
import subprocess
from collections import defaultdict
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


@TaskGenerator
def generate_yaml(sample_files, output_path, basedir):
    """Generate an NGLess-compatible YAML sample list."""
    lines = [f'basedir: {basedir}', 'samples:']
    for sample in sorted(sample_files):
        lines.append(f'  {sample}:')
        for fwd, rev in sample_files[sample]:
            lines.append('    - paired:')
            lines.append(f'        - {fwd}')
            lines.append(f'        - {rev}')

    with open(output_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return output_path


sample_files = defaultdict(list)

with open(BASE_DIR / 'RunInfo.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        sample = row['Title']
        urls = row['Download_path'].split('|')
        filenames = row['FileName'].split('|')
        pair = [None, None]
        for url, fname in zip(urls, filenames):
            opath = DATA_DIR / sample / fname
            download_file(url, str(opath))
            rel_path = str(Path(sample) / fname)
            if '_f1.' in fname:
                pair[0] = rel_path
            elif '_r2.' in fname:
                pair[1] = rel_path
        if pair[0] and pair[1]:
            sample_files[sample].append(tuple(pair))

generate_yaml(dict(sample_files), str(DATA_DIR / 'samples.yaml'), str(DATA_DIR))
