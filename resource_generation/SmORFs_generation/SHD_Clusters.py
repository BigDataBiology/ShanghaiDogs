## Map 100AA SmORFs to 90AA cluster representatives using CD-HIT results and export mappings to TSV
import os
import csv
import re
import fasta
from collections import defaultdict

#paths
base_dir = "/work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER"
input_fasta = os.path.join(base_dir, "100AA_SmORFs_sequences.faa")
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

# Map 100AA to 90AA IDs
id_mapping = {}
for cluster, seqs in cluster_mappings.items():
    representative_id = None
    for seq_id, similarity in seqs:
        if similarity == "*":
            representative_id = seq_id
            break
    if representative_id:
        rep_90aa_id = representative_id.replace("100AA", "90AA")
        for seq_id, similarity in seqs:
            id_mapping[seq_id] = (rep_90aa_id, similarity)

# Write output TSV
with open(output_tsv, "w", newline="") as out_f:
    writer = csv.writer(out_f, delimiter="\t")
    writer.writerow(["100AA SmORF ID", "90AA SmORF ID", "Similarity"])
    for seq_id_100aa in sorted(id_mapping.keys()):
        seq_id_90aa, similarity = id_mapping[seq_id_100aa]
        writer.writerow([seq_id_100aa, seq_id_90aa, similarity])

print(f"Output TSV created: {output_tsv}")
