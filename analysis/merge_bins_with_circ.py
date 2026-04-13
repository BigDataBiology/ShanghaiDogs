#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

base_binning = Path("/data/anna/animal_metagenome/long-mg-dog/04_binning/00_SemiBin2/LR")
base_assembly = Path("/data/anna/animal_metagenome/long-mg-dog/01_assembly")

all_samples = []
skipped_samples = []

for i in range(53):  # D000 to D052
    sample_id = f"D{i:03d}"

    bin_file = (
        base_binning
        / sample_id
        / "HQ_mq_MAGs"
        / "contig_bins_dastool.tsv"
    )

    assembly_file = (
        base_assembly
        / sample_id
        / "assembly_info.txt"
    )

    if not bin_file.exists():
        print(f"[WARNING] Missing bin file for {sample_id}, skipping.")
        skipped_samples.append(sample_id)
        continue

    if not assembly_file.exists():
        print(f"[WARNING] Missing assembly file for {sample_id}, skipping.")
        skipped_samples.append(sample_id)
        continue

    # ---- Read contig_bins_dastool.tsv ----
    bins_df = pd.read_csv(
        bin_file,
        sep="\t",
        header=None,
        names=["contig_id", "SemiBin_id"]
    )

    # REMOVE "_polypolish" SUFFIX
    bins_df["contig_id"] = bins_df["contig_id"].str.replace(
        r"_polypolish$", "", regex=True
    )

    bins_df["Sample_ID"] = sample_id

    # ---- Read assembly_info.txt ----
    assembly_df = pd.read_csv(
        assembly_file,
        sep="\t",
        comment="#",
        names=[
            "seq_name", "length", "cov", "circ",
            "repeat", "mult", "alt_group", "graph_path"
        ]
    )[["seq_name", "circ"]]

    # ---- Merge ----
    merged = bins_df.merge(
        assembly_df,
        left_on="contig_id",
        right_on="seq_name",
        how="left"
    ).drop(columns="seq_name")

    all_samples.append(merged)

# ---- Combine all samples ----
final_df = pd.concat(all_samples, ignore_index=True)

output_path = "/data/Projects/ShanghaiDogs/intermediate-outputs/tables/contig_bins_with_circ.tsv"
final_df.to_csv(output_path, sep="\t", index=False)

print(f"Final merged file written to: {output_path}")

if skipped_samples:
    print("!!! Skipped samples:")
    print(", ".join(skipped_samples))