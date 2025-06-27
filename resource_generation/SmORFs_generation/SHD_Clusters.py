### Clustering file

import os
import csv
import fasta
from tqdm import tqdm

# Paths
base_dir = "/work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER"
fasta_file = os.path.join(base_dir, "Mapped_SmORFs_sequences.faa")
cluster_file = os.path.join(base_dir, "clustered_SmORFs.clstr")
output_tsv = os.path.join(base_dir, "SHD_Clusters.tsv")
log_file = os.path.join(base_dir, "merge_log_Clusters.txt")

# Create dictionary from FASTA file (sequence -> 100AA ID)
seq_to_id_100aa = {}
for header, seq in tqdm(fasta.fasta_iter(fasta_file, full_header=True), desc="Processing FASTA file"):
    normalized_seq = seq.replace(" ", "").replace("\n", "").upper()
    seq_to_id_100aa[normalized_seq] = header  

# Create dictionary from cluster file (100AA ID -> 90AA representative ID)
id_100aa_to_90aa = {}
current_cluster_rep = None
with open(cluster_file, "r") as f:
    lines = f.readlines()
for line in tqdm(lines, desc="Processing cluster file"):
    line = line.strip()
    if line.startswith(">Cluster"):
        current_cluster_rep = None
    else:
        parts = line.split()
        smorf_id = parts[2].strip(">").split("...")[0]
        if parts[-1] == "*":
            current_cluster_rep = smorf_id  
        if current_cluster_rep:
            id_100aa_to_90aa[smorf_id.replace("90AA", "100AA")] = current_cluster_rep  

# Arrange sorted SmORF IDs by sequence
smorf_order = sorted(
    seq_to_id_100aa.items(),
    key=lambda x: x[0]  
)

with open(output_tsv, "w", newline='') as out_f:
    writer = csv.writer(out_f, delimiter="\t")
    writer.writerow(["100 AA SmORF ID", "90 AA SmORF ID"])
    for seq, smorf_id in tqdm(smorf_order, desc="Writing TSV"):
        writer.writerow([
            smorf_id,
            id_100aa_to_90aa.get(smorf_id, smorf_id.replace("100AA", "90AA"))
        ])

# Write log file
with open(log_file, "w") as log_f:
    log_f.write(f"Total unique sequences: {len(seq_to_id_100aa)}\n")

print(f"SHD Clusters TSV created: {output_tsv}")
print(f"Log file created: {log_file}")
