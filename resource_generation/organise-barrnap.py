import pandas as pd
from os import walk, makedirs
import gzip
meta = pd.read_csv('../data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv')
meta.set_index('Original ID', inplace=True)

TARGET_DIR = '../data/ShanghaiDogsMAGAnnotations/Barrnap/'
makedirs(TARGET_DIR, exist_ok=True)

found = set()
for base,_,fs in walk('../intermediate-outputs/07_ribosomal_genes/barrnap_fasta/'):
    for f in fs:
        if f.endswith('_ribosomal.fa'):
            f = f.removesuffix('_ribosomal.fa')
            sample, name = f.split('_', 1)
            fid = f"{name}_{sample}"
            if fid not in meta.index:
                continue
            if fid in found:
                print(f"Duplicate found for {fid}")
                raise SystemExit(1)
            found.add(fid)
            target_fname = meta.loc[fid, "Bin ID"].replace(".fna.gz", "_ribosomal.fna.gz")
            with gzip.open(f'{TARGET_DIR}/{target_fname}', 'wb') as out_f:
                with open(f'{base}/{f}_ribosomal.fa', 'rb') as in_f:
                    while chunk := in_f.read(32 * 1024):
                        out_f.write(chunk)

assert found == set(meta.index), 'Not all bins were found in the barrnap output!'
