#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import re

BASE = Path("/work/microbiome/shanghai_dogs")

SKANI_DIR = BASE / "resource_generation/MAGs_Onehealth/External_cohorts/Skani_lists/Quality_MAGs/Skani_Quality_Results"
SPIRE_META = BASE / "resource_generation/MAGs_Onehealth/External_cohorts/Skani_lists/spire_v1_genome_metadata.tsv"
SHD_META = BASE / "data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv"
HUMAN_ANI = BASE / "resource_generation/MAGs_Onehealth/SHD_Species_Rep_vs_Human_ani.tsv"

OUT_DIR = BASE / "resource_generation/MAGs_Onehealth"

EXT_OUT = OUT_DIR / "SHD_External_MAG_Matches.xlsx"
HUMAN_OUT = OUT_DIR / "SHD_CGMR_MAG_Matches.xlsx"

def clean_id(p):
    name = Path(str(p)).name
    return re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", name)

print("Loading metadata...")

spire = pd.read_csv(SPIRE_META, sep="\t", low_memory=False)
spire = spire[["genome_id", "classification"]].rename(columns={
    "genome_id": "spire_id",
    "classification": "spire_classification"
})

shd = pd.read_csv(SHD_META)
shd["shd_id"] = shd["Bin ID"].str.replace(".fna.gz", "", regex=False)
shd = shd[["shd_id", "Classification"]].rename(columns={
    "Classification": "shd_classification"
})

# EXTERNAL COHORT TABLE 
print("Processing external cohort files...")

dfs = []

for f in sorted(SKANI_DIR.glob("*_ani.tsv")):
    name = f.name.replace("_ani.tsv", "")
    parts = name.split("_")

    quality = parts[-1]
    if quality == "ALL":
        continue

    cohort = "_".join(parts[:-1])

    df = pd.read_csv(f, sep="\t")

    df["spire_id"] = df["Ref_file"].apply(clean_id)
    df["shd_id"] = df["Query_file"].apply(clean_id)
    df["cohort"] = cohort
    df["quality"] = quality

    df = df.rename(columns={
        "Align_fraction_ref": "cov_ref",
        "Align_fraction_query": "cov_query"
    })

    df = df[[
        "spire_id",
        "shd_id",
        "cohort",
        "quality",
        "ANI",
        "cov_ref",
        "cov_query"
    ]]

    dfs.append(df)

external_df = pd.concat(dfs, ignore_index=True)

external_df = external_df.merge(spire, on="spire_id", how="left")
external_df = external_df.merge(shd, on="shd_id", how="left")

external_df.to_excel(EXT_OUT, index=False)
print(f"Saved: {EXT_OUT}")

# HUMAN CGMR TABLE 
print("Processing human CGMR file...")

df = pd.read_csv(HUMAN_ANI, sep="\t")

df["gmr_id"] = df["Ref_file"].apply(clean_id)
df["shd_id"] = df["Query_file"].apply(clean_id)

df = df.rename(columns={
    "Align_fraction_ref": "cov_ref",
    "Align_fraction_query": "cov_query"
})

df = df[[
    "gmr_id",
    "shd_id",
    "ANI",
    "cov_ref",
    "cov_query"
]]

df = df.merge(shd, on="shd_id", how="left")

df.to_excel(HUMAN_OUT, index=False)

print(f"Saved: {HUMAN_OUT}")
