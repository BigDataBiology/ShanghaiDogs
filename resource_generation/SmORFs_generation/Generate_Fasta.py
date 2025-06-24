import os
from collections import defaultdict
import fasta  
import lib    

# Paths
base_dir = "/work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER"
output_fasta = os.path.join(base_dir, "Mapped_SmORFs_sequences.faa")
sequence_dict = defaultdict(list)
id_counter = 0

# Process each sample folder
for sample in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, sample)
    if os.path.isdir(sample_path):
        fasta_file = os.path.join(sample_path, "mapped.smorfs.faa")
        if os.path.exists(fasta_file):
            for header, seq in fasta.fasta_iter(fasta_file, full_header=True):
                sequence_dict[seq].append(header)

# Write unique sequences to output FASTA
with open(output_fasta, "w") as out_f:
    for seq, headers in sequence_dict.items():
        new_id = lib.pad9("SHD1_SM.100AA.", id_counter)  
        out_f.write(f">{new_id}\n{seq}\n")
        id_counter += 1 

print(f"FASTA file created: {output_fasta}")
