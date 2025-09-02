import os
import csv
from collections import defaultdict
import fasta
import lib
import gzip

# Paths
base_dir = "../../intermediate-outputs/GMSC_MAPPER"
habitat_taxonomy_tsv = os.path.join(base_dir, "100AA_SmORFs_habitat_taxonomy.tsv")
log_file = os.path.join(base_dir, "merge_log.txt")

os.makedirs(os.path.join(base_dir, "SHD_SMORF_resource"), exist_ok=True)
output_fasta = os.path.join(base_dir, "SHD_SMORF_resource", "100AA_SmORFs_sequences.faa.gz")
sequence_metadata_tsv = os.path.join(base_dir, "SHD_SMORF_resource", "100AA_SmORFs_origins.tsv.gz")


# all mapped SmORFs and group by sequence
sequence_dict = defaultdict(list)
for sample in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, sample)
    if os.path.isdir(sample_path):
        fasta_file = os.path.join(sample_path, "mapped.smorfs.faa")
        if os.path.exists(fasta_file):
            for header, seq in fasta.fasta_iter(fasta_file, full_header=True):
                clean_seq = seq.replace(" ", "").replace("\n", "").upper()
                sequence_dict[clean_seq].append(header)

# Assign IDs based on frequency (most frequent first)
sorted_sequences = sorted(sequence_dict.items(), key=lambda x: (-len(x[1]), x[0]))
seq_to_id = {}
with gzip.open(output_fasta, "wt") as out_f:
    for idx, (seq, headers) in enumerate(sorted_sequences):
        new_id = lib.pad6("SHD1_SM.100AA", idx)
        out_f.write(f">{new_id}\n{seq}\n")
        seq_to_id[seq] = new_id

print(f"FASTA file created: {output_fasta}")

#metadata, habitat, and taxonomy info
seq_metadata = defaultdict(list)
seq_counts = defaultdict(int)
smorf_to_habitat = defaultdict(set)
smorf_to_taxonomy = defaultdict(set)

for sample in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, sample)
    if not os.path.isdir(sample_path):
        continue
    sample_id = sample.replace("_PP1_PolcaCorr", "")

    pred_file = os.path.join(sample_path, "predicted.filterd.smorf.faa")
    if os.path.exists(pred_file):
        for header, seq in fasta.fasta_iter(pred_file, full_header=True):
            clean_seq = seq.replace(" ", "").replace("\n", "").upper()
            seq_counts[clean_seq] += 1

            # skip sequences not in GMSC (thus not in the mapped.faa file)
            if clean_seq not in seq_to_id:
                continue

            parts = header.split("#")
            assert len(parts) >= 5

            smorf_id = parts[0].strip()
            contig = parts[1].strip()
            start = parts[2].strip()
            end = parts[3].strip()
            strand = "+" if parts[4].strip() == "1" else "-"
            coord = f"{start}-{end}"
            seq_metadata[clean_seq].append({
                "sample_id": sample_id,
                "contig": contig,
                "coords": coord,
                "strand": strand,
                "smorf_id": smorf_id
            })

    # habitat.out.smorfs.tsv
    habitat_file = os.path.join(sample_path, "habitat.out.smorfs.tsv")

    if os.path.exists(habitat_file):
        with open(habitat_file, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader, None)
            for row in reader:
                smorf_id = row[0].strip()
                habitats = [h.strip() for h in row[1].split(",") if h.strip()]
                assert (smorf_id, sample_id) not in smorf_to_habitat
                smorf_to_habitat[(smorf_id, sample_id)] = set(habitats)

    #taxonomy.out.smorfs.tsv
    taxonomy_file = os.path.join(sample_path, "taxonomy.out.smorfs.tsv")
    if os.path.exists(taxonomy_file):
        with open(taxonomy_file, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader, None)
            for row in reader:
                smorf_id = row[0].strip()
                taxonomy = row[1].strip()
                assert (smorf_id, sample_id) not in smorf_to_taxonomy
                smorf_to_taxonomy[(smorf_id, sample_id)] = taxonomy

#Order sequences by frequency and write metadata
smorf_order = sorted(seq_to_id.items(), key=lambda x: (-seq_counts[x[0]], x[0]))

# SmORFs_metadata.tsv
with gzip.open(sequence_metadata_tsv, "wt", newline='') as out_f:
    writer = csv.writer(out_f, delimiter="\t")
    writer.writerow(["SmORF ID", "Sample ID", "Contig", "Coordinates", "Strand"])
    for seq, smorf_id in smorf_order:
        for entry in seq_metadata[seq]:
            clean_contig = entry["contig"].rsplit("_", maxsplit=1)[0]
            writer.writerow([
                smorf_id,
                entry["sample_id"],
                clean_contig,
                entry["coords"],
                entry["strand"]
            ])

# SmORFs_habitat_taxonomy.tsv
with open(habitat_taxonomy_tsv, "w", newline='') as out_f:
    writer = csv.writer(out_f, delimiter="\t")
    writer.writerow(["SmORF ID", "Habitat", "Taxonomy"])
    for seq, smorf_id in smorf_order:
        seen_samples = set()
        habitats = []
        taxonomies = []
        for entry in seq_metadata[seq]:
            sample_id = entry["sample_id"]
            smorf_id = entry["smorf_id"]
            habitats.append(smorf_to_habitat[(smorf_id, sample_id)])
            taxonomies.append(smorf_to_taxonomy[(smorf_id, sample_id)])
        habitat = habitats[0]
        taxonomy = taxonomies[0]
        assert all(habitat == habitats[0] for h in habitats), "Inconsistent habitats for the same smORF"
        assert all(taxonomy == taxonomies[0] for t in taxonomies), "Inconsistent taxonomies for the same smORF"
        habitat = ','.join(sorted(habitat))
        writer.writerow([smorf_id, sample_id, habitat, taxonomy])

# Step 5: Log file
with open(log_file, "w") as log_f:
    log_f.write(f"Total unique sequences: {len(seq_to_id)}\n")
    duplicate_count = sum(1 for count in seq_counts.values() if count > 1)
    log_f.write(f"Number of sequences with duplicates: {duplicate_count}\n")
    log_f.write("Detailed duplicate counts with associated smORFs:\n")
    for seq, smorf_id in smorf_order:
        count = seq_counts[seq]
        if count > 1:
            smorfs = ", ".join(sorted(set(entry["smorf_id"] for entry in seq_metadata[seq])))
            log_f.write(f"Sequence with ID {smorf_id} appeared {count} times (smORFs: {smorfs})\n")

print(f" FASTA created: {output_fasta}")
print(f" Metadata TSV: {sequence_metadata_tsv}")
print(f" Habitat & Taxonomy TSV: {habitat_taxonomy_tsv}")
print(f" Log: {log_file}")
