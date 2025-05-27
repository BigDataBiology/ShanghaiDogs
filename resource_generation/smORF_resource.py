#Running GMSC MAPPER on All SHD dogs Assemblies data 
conda activate gmscmapper

# paths
INPUT_DIR='shanghai_dogs\data\ShanghaiDogsAssemblies'
OUTPUT_DIR='shanghai_dogs/intermediate-outputs/GMSC_MAPPER_TEST'
GMSC='shanghai_dogs/intermediate-outputs/GMSC_Mapper_dogs/GMSC-mapper'

# Loop over all files with .fa, .fasta, or .fna extensions in the input directory
for file in "$INPUT_DIR"/*.{fa,fasta,fna}; do
    # Skip if no files are found
    [ -e "$file" ] || continue

    filename=$(basename "$file")
    filename_no_ext="${filename%.*}"
    outdir="$OUTPUT_DIR/$filename_no_ext"

    if [ -d "$outdir" ] && [ "$(ls -A "$outdir")" ]; then
        echo "✅ Skipping $filename (output exists)"
        continue
    fi

    echo "🚀 Running gmsc-mapper on $filename"
    mkdir -p "$outdir"
    gmsc-mapper -i "$file" -o "$outdir" --dbdir "$GMSC/db" -t 24
done

### PART 2 Generates a TSV File
import pandas as pd
import os
import gzip
import re

# Define paths
base_dir = "/work/microbiome/shanghai_dogs"
sample_dir = "D030_PP1_PolcaCorr"
sample_number = re.match(r"D(\d+)_PP1_PolcaCorr", sample_dir).group(1)

smorf_file = f"{base_dir}/intermediate-outputs/GMSC_MAPPER_TEST/{sample_dir}/predicted.filterd.smorf.faa"
alignment_file = f"{base_dir}/intermediate-outputs/GMSC_MAPPER_TEST/{sample_dir}/alignment.out.smorfs.tsv"
habitat_file = f"{base_dir}/intermediate-outputs/GMSC_MAPPER_TEST/{sample_dir}/habitat.out.smorfs.tsv"
mag_dir = f"{base_dir}/data/ShanghaiDogsMAGs"
mimag_report = f"{base_dir}/data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv"
output_tsv = f"{base_dir}/intermediate-outputs/GMSC_MAPPER_TEST/{sample_dir}/D030_smorfs.tsv"

# Step 1: Parse smORF IDs and contig IDs from predicted.filterd.smorf.faa
smorfs = []
smorf_count = 1
with open(smorf_file, "r") as f:
    for line in f:
        if line.startswith(">"):
            parts = line[1:].strip().split("#")
            orig_id = parts[0].strip()
            contig_id = parts[1].strip()
            core_contig = "_".join(contig_id.split("_")[:-1])
            new_id = f"SHD_SM.{sample_number}.{smorf_count:05d}"
            smorfs.append({
                "orig_id": orig_id,
                "new_id": new_id,
                "contig_id": contig_id,
                "core_contig": core_contig
            })
            smorf_count += 1

smorfs_df = pd.DataFrame(smorfs)

# Step 2: Parse GMSC hits
gmsc_hits = {}
if os.path.exists(alignment_file):
    alignment_df = pd.read_csv(alignment_file, sep="\t", header=None, usecols=[0, 1], names=["qseqid", "sseqid"])
    for qseqid, group in alignment_df.groupby("qseqid"):
        hits = ";".join(group["sseqid"])
        gmsc_hits[qseqid] = hits
smorfs_df["gmsc_hits"] = smorfs_df["orig_id"].map(gmsc_hits).fillna("-")

# Step 3: Parse habitats
habitats = {}
if os.path.exists(habitat_file):
    habitat_df = pd.read_csv(habitat_file, sep="\t", names=["qseqid", "habitat"])
    for _, row in habitat_df.iterrows():
        habitats[row["qseqid"]] = row["habitat"]
smorfs_df["habitat"] = smorfs_df["orig_id"].map(habitats).fillna("-")

# Step 4: Map contigs to MAG IDs
mag_mapping = {}
mag_files = [f for f in os.listdir(mag_dir) if f.endswith(".fna.gz")]
for mag_file in mag_files:
    with gzip.open(os.path.join(mag_dir, mag_file), "rt") as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                if header.endswith(" D030"):
                    parts = header.split()
                    mag_id = parts[0].split("_")[0] + "_" + parts[0].split("_")[1]
                    contig_id = parts[1]
                    core_contig = contig_id
                    mag_mapping[core_contig] = mag_id

smorfs_df["mag_id"] = smorfs_df["core_contig"].map(mag_mapping).fillna("-")
unmatched = smorfs_df[smorfs_df["mag_id"] == "-"]["core_contig"].unique()
for contig in unmatched:
    print(f"Unmatched smORF contig: {contig}")

# Step 5: Map MAG ID to taxonomy
mimag_df = pd.read_csv(mimag_report)
taxonomy_mapping = {}
for _, row in mimag_df.iterrows():
    bin_id = row["Bin ID"]
    mag_id = bin_id.replace(".fna.gz", "")
    taxonomy_mapping[mag_id] = row["Classification"]

smorfs_df["taxonomy"] = smorfs_df["mag_id"].map(taxonomy_mapping).fillna("-")

# Step 6: Write TSV with renamed IDs
output_df = pd.DataFrame({
    "ID": smorfs_df["new_id"],
    "GMSC_Hits": smorfs_df["gmsc_hits"],
    "Taxonomy": smorfs_df["taxonomy"],
    "Other_Habitats": smorfs_df["habitat"]
})

output_df.to_csv(output_tsv, sep="\t", index=False)
print(f"TSV saved to {output_tsv}")
print(f"Total smORFs processed: {len(smorfs_df)}")

### Part3 Generates a Fasta File
#Paths

base_dir = "/work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER_TEST"
sample_dir = "D030_PP1_PolcaCorr"
input_faa = f"{base_dir}/{sample_dir}/mapped.smorfs.faa"
output_faa = f"{base_dir}/{sample_dir}/SHD_mapped_smorfs.faa"

# Extract
sample_number = re.match(r"D(\d+)_PP1_PolcaCorr", sample_dir).group(1)

# Renaming smORFs
with open(input_faa, "r") as f_in, open(output_faa, "w") as f_out:
    smorf_count = 1
    for line in f_in:
        if line.startswith(">"):
            # Original ID: >smORF_00051
            new_id = f">SHD_SM.{sample_number}.{smorf_count:05d}"
            f_out.write(new_id + "\n")
            smorf_count += 1
        else:
            # Write sequence unchanged
            f_out.write(line)

print(f"Renamed smORFs saved to {output_faa}")
print(f"Total smORFs renamed: {smorf_count - 1}")
