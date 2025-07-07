## Map 100AA SmORFs to 90AA cluster representatives using CD-HIT results and export mappings to TSV
import os
import gzip
import re
import fasta
from collections import defaultdict
import lib
import pandas as pd

#paths
base_dir = "/work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER/"
input_fasta = os.path.join(base_dir, "100AA_SmORFs_sequences.faa.gz")
origin_tsv = os.path.join(base_dir, "100AA_SmORFs_origins.tsv.gz")
cluster_file = os.path.join(base_dir, "clustered_SmORFs.clstr")
output_tsv = os.path.join(base_dir, "SHD_Clusters.tsv")

# Parse CD-HIT cluster file
def parse_cdhit_clusters(cluster_file):
    cluster_mappings = defaultdict(list)
    current_cluster = None
    with open(cluster_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                current_cluster = line
            else:
                match = re.match(r"\d+\s+\d+aa,\s+>(.*?)\.\.\.\s*(.*)", line)
                if match:
                    seq_id = match.group(1)
                    similarity = match.group(2).strip()
                    if similarity == "*":
                        cluster_mappings[current_cluster].append((seq_id, "*"))
                    else:
                        perc_match = re.search(r"at\s+(\d+\.\d+%)", similarity)
                        if perc_match:
                            percentage = perc_match.group(1)
                            cluster_mappings[current_cluster].append((seq_id, percentage))
    return cluster_mappings

#FASTA file
fasta_sequences = {header: seq for header, seq in fasta.fasta_iter(input_fasta, full_header=True)}

# cluster file
cluster_mappings = parse_cdhit_clusters(cluster_file)

# Verify sequence IDs
cluster_ids = set(seq_id for seqs in cluster_mappings.values() for seq_id, _ in seqs)
fasta_ids = set(fasta_sequences.keys())
missing_ids = cluster_ids - fasta_ids
if missing_ids:
    raise ValueError(f"Sequence IDs in .clstr not found in FASTA: {missing_ids}")

origins = pd.read_csv(origin_tsv, sep="\t")
count100aa = origins['SmORF ID'].value_counts().to_dict()

def count_smorfs(mapping_item):
    # Sort by the number of SmORFs in the cluster and then by the minimum 100AA
    # SmORF ID as a tiebreaker
    return (
            sum(
                count100aa[seq_id] for seq_id, _ in mapping_item[1]
                ),
            min(seq_id for seq_id, _ in mapping_item[1])
            )

clusters = list(cluster_mappings.items())
clusters.sort(key=count_smorfs, reverse=True)
# Map 100AA to 90AA IDs
mappings = []
for idx, (_, seqs) in enumerate(clusters):
    rep_id = lib.pad6("SHD1_SM.90AA", idx)
    for seq_id, similarity in seqs:
        mappings.append((seq_id, rep_id, similarity))

mappings = pd.DataFrame(mappings, columns=["100AA SmORF ID", "90AA SmORF ID", "Similarity"])
mappings.sort_values(by=["100AA SmORF ID"], inplace=True)
mappings.to_csv(output_tsv, sep="\t", index=False)

assert mappings.eval('Similarity == "*"').sum() == len(set(mappings['90AA SmORF ID']))
mappings = mappings.query('Similarity == "*"')
assert len(set(mappings['90AA SmORF ID'])) == len(set(mappings['100AA SmORF ID']))
map_90aa_100aa = mappings.set_index('100AA SmORF ID')['90AA SmORF ID'].to_dict()

seqs = []
for h,seq in fasta.fasta_iter(input_fasta):
    if h in map_90aa_100aa:
        seqs.append((map_90aa_100aa[h], h, seq))
assert len(seqs) == len(map_90aa_100aa)

seqs.sort()

with open(f"{base_dir}/SHD_90AA_SmORFs.faa", "wt") as out_f:
    for h90, h100, seq in seqs:
        out_f.write(f">{h90} {h100}\n{seq}\n")

