from jug import TaskGenerator, bvalue
import requests
from glob import glob
import fasta
from collections import defaultdict


@TaskGenerator
def lookup_1seq(h, seq):
    from time import sleep
    sleep(0.2)  # to avoid hitting the server too hard
    URL = 'https://microbeatlas.org/ma_gw/index.py'
    payload = {
        "action": "mapseq",
        "seq": f'>{h}\\n{seq}',
        "seq_original": f'>{h}\n{seq}',
        "command": f'@map:mapseq(">{h}\\n{seq}", "-1")',
    }
    payload

    r = requests.post(
            URL,
            data=payload,
            headers = {
                'User-Agent': 'Python/requests (custom script. Hi from Brisbane, this is Luis. How is everyone doing in Zurich?)',
            }
        )
    if r.status_code != 200:
        print(r.status_code)
        raise OSError(f"Error {r.status_code} from Microbe Atlas server")
    return r.status_code, r.json()

@TaskGenerator
def get_sequences():
    """
    Get sequences from the Shanghai Dogs MAG Annotations.
    """
    fs = sorted(glob('../data/ShanghaiDogsMAGAnnotations/Barrnap/*.fna.gz'))
    seqs = defaultdict(list)
    for f in fs:
        binid = f.split('/')[-1].removesuffix('_ribosomal.fna.gz')  # just to get the filename
        for h,seq in fasta.fasta_iter(f):
            if h.startswith('16S_rRNA::'):
                h = h.split('::')[1]
                seqs[seq].append((binid, h))
    return seqs

results = []
for seq, hs in bvalue(get_sequences()).items():
    h = hs[0]
    results.append((seq, hs, lookup_1seq(h, seq)))
