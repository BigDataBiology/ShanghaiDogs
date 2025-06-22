import pandas as pd
import os

#paths
mimag_report_path = 'data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv'
args_filt_path = 'intermediate-outputs/06_ARG/MAGs-ARGs_ALL_filt.txt'

mimag_df = pd.read_csv(mimag_report_path)
args_df = pd.read_csv(args_filt_path)

# Extract species
mimag_df['Species'] = mimag_df['Classification'].str.extract(r's__([^;]+)')
base_df = mimag_df[['Bin ID', 'Species']].copy()

# potential pathogens to filter
pathogens = [
    'Campylobacter_D upsaliensis', 'Enterobacter hormaechei_A', 'Enterobacter roggenkampii',
    'Enterococcus faecalis', 'Enterococcus_B faecium', 'Enterococcus_B hirae',
    'Enterococcus_B lactis', 'Enterococcus_D gallinarum', 'Helicobacter_A bilis_A',
    'Helicobacter_A rappini', 'Helicobacter_B canis', 'Helicobacter_B canis_B',
    'Helicobacter_C magdeburgensis', 'Klebsiella pneumoniae', 'Proteus mirabilis'
]

#base_df to include only potential pathogens
base_df = base_df[base_df['Species'].isin(pathogens)]

# Count ARGs per Bin
arg_counts = args_df.pivot_table(
    index='Bin ID',
    columns='Best_Hit_ARO',
    aggfunc='size',
    fill_value=0
).reset_index()

full_df = pd.merge(base_df, arg_counts, on='Bin ID', how='left')

# Fill bins with no ARGs with 0
full_df = full_df.fillna(0)

for col in full_df.columns[2:]:
    full_df[col] = full_df[col].astype(int)

# Collapse by species, calculating median for each ARG
collapsed_df = full_df.groupby('Species').median(numeric_only=True).reset_index()

# Add total ARGs column - sum of medians across ARGs for each species
collapsed_df['Total_ARGs'] = collapsed_df.iloc[:, 1:].sum(axis=1)

print(collapsed_df)
