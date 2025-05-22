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

###PART 2
import pandas as pd
import os
import gzip

base_dir = "/work/microbiome/shanghai_dogs"
smorf_file = f"{base_dir}/intermediate-outputs/GMSC_MAPPER_TEST/D030_PP1_PolcaCorr/predicted.filterd.smorf.faa"
habitat_file = f"{base_dir}/intermediate-outputs/GMSC_MAPPER_TEST/D030_PP1_PolcaCorr/habitat.out.smorfs.tsv"
mag_dir = f"{base_dir}/data/ShanghaiDogsMAGs"
mimag_report = f"{base_dir}/data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv"
output_tsv = f"{base_dir}/intermediate-outputs/GMSC_MAPPER_TEST/D030_PP1_PolcaCorr/D030_smorfs.tsv"

# Step 1: Parse smORF IDs and contig IDs from predicted.filtered.smorf.faa
smorfs = []
with open(smorf_file, "r") as f:
    for line in f:
        if line.startswith(">"):
            parts = line[1:].strip().split("#")
            smorf_id = parts[0].strip()  # e.g., smORF_00003
            contig_id = parts[1].strip()  # e.g., contig_1179_polypolish_5
            core_contig = "_".join(contig_id.split("_")[:-1])  # e.g., contig_1179_polypolish
            smorfs.append({"smorf_id": smorf_id, "contig_id": contig_id, "core_contig": core_contig})

smorfs_df = pd.DataFrame(smorfs)

# Step 2: Parse habitats from habitat.out.smorfs.tsv
habitats = {}
if os.path.exists(habitat_file):
    habitat_df = pd.read_csv(habitat_file, sep="\t", names=["qseqid", "habitat"])
    for _, row in habitat_df.iterrows():
        habitats[row["qseqid"]] = row["habitat"]
smorfs_df["habitat"] = smorfs_df["smorf_id"].map(habitats).fillna("-")

# Step 3: Map contig IDs to MAG IDs
mag_mapping = {}
mag_files = [f for f in os.listdir(mag_dir) if f.endswith(".fna.gz")]
for mag_file in mag_files:
    with gzip.open(os.path.join(mag_dir, mag_file), "rt") as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                if " D030" in header:
                    parts = header.split()
                    mag_id = parts[0].split("_")[0] + "_" + parts[0].split("_")[1]  # e.g., SHD1_1440
                    contig_id = parts[1] 
                    core_contig = contig_id  
                    mag_mapping[core_contig] = mag_id


smorfs_df["mag_id"] = smorfs_df["core_contig"].map(mag_mapping).fillna("-")

# Step 4: Load MIMAG report and map MAG IDs to taxonomy
mimag_df = pd.read_csv(mimag_report)
taxonomy_mapping = {}
for _, row in mimag_df.iterrows():
    bin_id = row["Bin ID"]
    mag_id = bin_id.replace(".fna.gz", "")
    taxonomy_mapping[mag_id] = row["Classification"]

smorfs_df["taxonomy"] = smorfs_df["mag_id"].map(taxonomy_mapping).fillna("-")

# Step 5: Create TSV with required columns
output_df = pd.DataFrame({
    "ID": smorfs_df["smorf_id"],
    "GMSC_Hits": "-",
    "Taxonomy": smorfs_df["taxonomy"],
    "Other_Habitats": smorfs_df["habitat"]
})

output_df.to_csv(output_tsv, sep="\t", index=False)
print(f"TSV saved to {output_tsv}")
