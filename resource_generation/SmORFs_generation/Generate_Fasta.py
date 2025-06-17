import os
from collections import defaultdict

# paths
base_dir = "/work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER"
output_fasta = os.path.join(base_dir, "Mapped_SmORFs_sequences.faa")
sequence_dict = defaultdict(list)
id_counter = 0

# parse FASTA file (I will add the symlink soon had issue with auhtnetication)
def fasta_iter(fname, full_header=False):
    header = None
    chunks = []
    if fname.endswith('.gz'):
        import gzip
        op = gzip.open
    elif fname.endswith('.bz2'):
        import bz2
        op = bz2.open
    elif fname.endswith('.xz'):
        import lzma
        op = lzma.open
    else:
        op = open
    with op(fname, 'rt') as f:
        for line in f:
            if line[0] == '>':
                if header is not None:
                    yield header, ''.join(chunks)
                line = line[1:].strip()
                if not line:
                    header = ''
                elif full_header:
                    header = line.strip()
                else:
                    header = line.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, ''.join(chunks)

# Process each sample folder
for sample in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, sample)
    if os.path.isdir(sample_path):
        fasta_file = os.path.join(sample_path, "mapped.smorfs.faa")
        if os.path.exists(fasta_file):
            for header, seq in fasta_iter(fasta_file, full_header=True):
                sequence_dict[seq].append(header)

# unique sequences to output FASTA
with open(output_fasta, "w") as out_f:
    for seq, headers in sequence_dict.items():
        new_id = f"SHD.SM.100AA.000_{id_counter:03d}"
        out_f.write(f">{new_id}\n{seq}\n")
        id_counter += 1

print(f"FASTA file created: {output_fasta}")
