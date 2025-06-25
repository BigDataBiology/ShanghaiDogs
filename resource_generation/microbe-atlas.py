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


@TaskGenerator
def reorganize_results(results):
    """
    Reorganize the results from the Microbe Atlas server into a DataFrame.
    Parameters:
    results (list): List of tuples containing sequence, headers, and response from the server.

    Returns:
    pd.DataFrame: DataFrame with columns ['MAG', 'Seq', 'OTU'].
    """
    import pandas as pd
    reorg = []
    for seq, hs, (st, r) in results:
        assert st == 200, f"Error {st} from Microbe Atlas server"
        otus = r['otutable']['otus']
        if len(otus) == 0:
            otus = 'NA'
        elif len(otus) == 1:
            [otu] = otus
        else:
            raise ValueError(f"Unexpected number of OTUs: {len(otus)} for sequence {seq}")
        for m, _ in hs:
            reorg.append((m, seq, otu))
    data = pd.DataFrame(reorg, columns=['MAG', 'Seq', 'OTU'])
    data.sort_values(by='MAG', inplace=True)
    return data.reset_index(drop=True)

@TaskGenerator
def save_results(final):
    """
    Save the final DataFrame to a CSV file.
    """
    import pandas as pd
    final.to_csv('../intermediate-outputs/07_ribosomal_genes/microbe-atlas.csv.gz', index=False)

results = []
for seq, hs in bvalue(get_sequences()).items():
    h = hs[0][1]
    results.append((seq, hs, lookup_1seq(h, seq)))
final = reorganize_results(results)
save_results(final)
