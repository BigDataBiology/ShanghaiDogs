import os
import pandas as pd
import fasta

# Paths
base_dir = "/work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER"
input_fasta = os.path.join(base_dir, "100AA_SmORFs_sequences.faa")
cluster_tsv = os.path.join(base_dir, "SHD_Clusters.tsv")

clusters = pd.read_csv(cluster_tsv, sep="\t")

assert clusters.eval('Similarity == "*"').sum() == len(set(clusters['90AA SmORF ID']))
clusters = clusters.query('Similarity == "*"')
assert len(set(clusters['90AA SmORF ID'])) == len(set(clusters['100AA SmORF ID']))
map_90aa_100aa = clusters.set_index('100AA SmORF ID')['90AA SmORF ID'].to_dict()

seqs = []
for h,seq in fasta.fasta_iter(input_fasta):
    if h in map_90aa_100aa:
        seqs.append((map_90aa_100aa[h], h, seq))
assert len(seqs) == len(map_90aa_100aa)

seqs.sort()

(f"{base_dir}/SHD_90AA_SmORFs.faa", "wt")
with open(f"{base_dir}/SHD_90AA_SmORFs.faa", "wt") as out_f:
    for h90, h100, seq in seqs:
        out_f.write(f">{h90} {h100}\n{seq}\n")

