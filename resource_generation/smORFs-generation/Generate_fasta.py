## FASTA FILE FOR SHANGHAI DOG SMORFS
import os
import re

# GMSC-Mapper output
GMSC_output_data = "/work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER"
sample_dirs = sorted([d for d in os.listdir(GMSC_output_data) if re.match(r"D\d{3}_PP1_PolcaCorr", d)])

for sample_dir in sample_dirs:
    # Input and output FASTA paths
    input_faa = os.path.join(GMSC_output_data, sample_dir, "mapped.smorfs.faa")
    output_faa = os.path.join(GMSC_output_data, sample_dir, "SHD_mapped_smorfs.faa")

    # Process FASTA files to :SHD_SM.100AA.<count>
    with open(input_faa, "r") as f_in, open(output_faa, "w") as f_out:
        smorf_count = 1
        for line in f_in:
            if line.startswith(">"):
                f_out.write(f">SHD_SM.100AA.{smorf_count:05d}\n")
                smorf_count += 1
            else:
                f_out.write(line)

    print(f"[DONE] {sample_dir}: Renamed {smorf_count - 1} smORFs → {output_faa}")

## TSV FILEs For SHANGHAI DOG SMORFS

import os
import re
import pandas as pd

# GMSC-Mapper output directory
GMSC_output_data = "/work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER"
sample_dirs = sorted([d for d in os.listdir(GMSC_output_data) if re.match(r"D\d{3}_PP1_PolcaCorr", d)])

for sample_dir in sample_dirs:
    # Define file paths
    habitat_file = os.path.join(GMSC_output_data, sample_dir, "habitat.out.smorfs.tsv")
    taxonomy_file = os.path.join(GMSC_output_data, sample_dir, "taxonomy.out.smorfs.tsv")
    faa_file = os.path.join(GMSC_output_data, sample_dir, "predicted.filterd.smorf.faa")
    info_tsv = os.path.join(GMSC_output_data, sample_dir, "SHD_smorfs_info.tsv")
    annotation_tsv = os.path.join(GMSC_output_data, sample_dir, "SHD_smorfs.tsv")

    habitat_df = pd.read_csv(habitat_file, sep="\t")
    taxonomy_df = pd.read_csv(taxonomy_file, sep="\t")
    habitat_df = habitat_df.rename(columns={"qseqid": "smORF_ID", "habitat": "Habitat"})
    taxonomy_df = taxonomy_df.rename(columns={"q_seqid": "smORF_ID", "taxonomy": "Taxonomy"})

    # Parse FASTA file 
    faa_data = []
    with open(faa_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                # Extract smORF ID, contig, coordinates, and strand from header
                match = re.match(r">(\S+) # (\S+) # (\d+) # (\d+) # ([+-]1)", line)
                if match:
                    smorf_id = match.group(1)
                    contig = match.group(2)
                    start = match.group(3)
                    end = match.group(4)
                    strand = match.group(5)
                    smorf_number = re.search(r"\d+$", smorf_id).group()
                    # Create new smORF ID with SHD format using smORF_ID number
                    new_smorf_id = f"SHD_SM.100AA.{smorf_number.zfill(5)}"
                    faa_data.append({
                        "smORF_ID": smorf_id,
                        "New_smORF_ID": new_smorf_id,
                        "Contig": contig,
                        "Start": start,
                        "End": end,
                        "Strand": strand
                    })

    faa_df = pd.DataFrame(faa_data)

    # Create info DataFrame (original smORF_ID, new smORF_ID, Contig, Start, End, Strand)
    info_df = faa_df[["smORF_ID", "New_smORF_ID", "Contig", "Start", "End", "Strand"]]
    info_df = info_df.rename(columns={"New_smORF_ID": "Sample_ID"})

    # Merge DataFrames for annotation
    merged_df = faa_df.merge(habitat_df, on="smORF_ID", how="left")
    merged_df = merged_df.merge(taxonomy_df, on="smORF_ID", how="left")

    # Create annotation DataFrame
    annotation_df = merged_df[["New_smORF_ID", "Habitat", "Taxonomy"]]
    annotation_df = annotation_df.rename(columns={"New_smORF_ID": "Sample_ID"})

    # Save TSV files
    info_df.to_csv(info_tsv, sep="\t", index=False)
    annotation_df.to_csv(annotation_tsv, sep="\t", index=False)

    print(f"[DONE] {sample_dir}: Generated {info_tsv} with {len(info_df)} smORFs")
    print(f"[DONE] {sample_dir}: Generated {annotation_tsv} with {len(annotation_df)} smORFs")
