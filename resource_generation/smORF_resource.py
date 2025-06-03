import pandas as pd
import os
import gzip
import re

# Define paths relative to working directory
mag_dir = "data/ShanghaiDogsMAGs"
mimag_report = "data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv"
gmsc_output_dir = "intermediate-outputs/GMSC_MAPPER_TEST"

# Get all sample directories
sample_dirs = sorted([d for d in os.listdir(gmsc_output_dir) if re.match(r"D\d{3}_PP1_PolcaCorr", d)])

for sample_dir in sample_dirs:
    # Extract sample number
    sample_number = re.match(r"D(\d+)_PP1_PolcaCorr", sample_dir).group(1)

    # Define file paths
    smorf_file = f"{gmsc_output_dir}/{sample_dir}/predicted.filterd.smorf.faa"
    alignment_file = f"{gmsc_output_dir}/{sample_dir}/alignment.out.smorfs.tsv"
    habitat_file = f"{gmsc_output_dir}/{sample_dir}/habitat.out.smorfs.tsv"
    output_tsv = f"{gmsc_output_dir}/{sample_dir}/SHD_D{sample_number}_smorfs.tsv"

    # Step 1: Parse smORF IDs and contig IDs
    smorfs = []
    smorf_count = 1
    with open(smorf_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                parts = line[1:].strip().split("#")
                if len(parts) < 2:
                    raise ValueError(f"Invalid header format in {smorf_file}: {line}")
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
    alignment_df = pd.read_csv(alignment_file, sep="\t", header=None, usecols=[0, 1], names=["qseqid", "sseqid"])
    gmsc_hits = alignment_df.groupby("qseqid")["sseqid"].apply(lambda x: ";".join(x)).to_dict()
    smorfs_df["gmsc_hits"] = smorfs_df["orig_id"].map(gmsc_hits).fillna("-")

    # Step 3: Parse habitats
    habitat_df = pd.read_csv(habitat_file, sep="\t", names=["qseqid", "habitat"])
    habitats = dict(zip(habitat_df["qseqid"], habitat_df["habitat"]))
    smorfs_df["habitat"] = smorfs_df["orig_id"].map(habitats).fillna("-")

    # Step 4: Map contigs to MAG IDs
    mag_mapping = {}
    mag_files = [f for f in os.listdir(mag_dir) if f.endswith(".fna.gz")]
    for mag_file in mag_files:
        with gzip.open(os.path.join(mag_dir, mag_file), "rt") as f:
            for line in f:
                if line.startswith(">"):
                    header = line[1:].strip()
                    if f"D{sample_number}" in header:
                        parts = header.split()
                        mag_id = parts[0].split("_")[0] + "_" + parts[0].split("_")[1]
                        contig_id = parts[1]
                        core_contig = contig_id
                        mag_mapping[core_contig] = mag_id

    smorfs_df["mag_id"] = smorfs_df["core_contig"].map(mag_mapping).fillna("-")
    unmatched = smorfs_df[smorfs_df["mag_id"] == "-"]["core_contig"].unique()
    for contig in unmatched:
        print(f"[INFO] {sample_dir}: Unmatched smORF contig: {contig}")

    # Step 5: Map MAG ID to taxonomy
    mimag_df = pd.read_csv(mimag_report)
    taxonomy_mapping = dict(zip(
        mimag_df["Bin ID"].str.replace(".fna.gz", "", regex=False),
        mimag_df["Classification"]
    ))
    smorfs_df["taxonomy"] = smorfs_df["mag_id"].map(taxonomy_mapping).fillna("-")

    # Step 6: Write TSV with renamed IDs
    output_df = pd.DataFrame({
        "ID": smorfs_df["new_id"],
        "GMSC_Hits": smorfs_df["gmsc_hits"],
        "Taxonomy": smorfs_df["taxonomy"],
        "Other_Habitats": smorfs_df["habitat"]
    })

    output_df.to_csv(output_tsv, sep="\t", index=False)
    print(f"[DONE] {sample_dir}: TSV saved to {output_tsv}")
    print(f"[DONE] {sample_dir}: Total smORFs processed: {len(smorfs_df)}")


##GMSC OUTPUT
#Running GMSC MAPPER on All SHD dogs Assemblies data 
conda activate gmscmapper

# paths
INPUT_DIR='shanghai_dogs/data/ShanghaiDogsAssemblies'
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

