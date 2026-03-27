import pandas as pd
from os import walk, makedirs
import gzip
import fasta
from pathlib import Path
meta = pd.read_csv('../data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv')
meta.set_index('Original ID', inplace=True)

TARGET_DIR = '../data/ShanghaiDogsMAGAnnotations/Barrnap/'
FNA_DIR = Path('../data/ShanghaiDogsMAGs/')
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
            fna_name = FNA_DIR / f"{meta.loc[fid, 'Bin ID']}"
            renames = {}
            for h,_ in fasta.fasta_iter(str(fna_name), full_header=True):
                new,old,_ = h.split(' ', 2)
                renames[old] = new

            with gzip.open(f'{TARGET_DIR}/{target_fname}', 'wb') as out_f:
                for h, seq in fasta.fasta_iter(f'{base}/{f}_ribosomal.fa', full_header=True):
                    contig = h.split(':')[2]
                    new_h = h.replace(contig, renames[contig])
                    out_f.write(f'>{new_h}\n'.encode())
                    out_f.write(f'{seq}\n'.encode())

assert found == set(meta.index), 'Not all bins were found in the barrnap output!'
