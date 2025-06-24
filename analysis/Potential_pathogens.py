import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
mimag_report_path = 'data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv'
args_filt_path = 'intermediate-outputs/06_ARG/MAGs-ARGs_ALL_filt.txt'

# Read data
mimag_df = pd.read_csv(mimag_report_path)
args_df = pd.read_csv(args_filt_path)

# Extract species
mimag_df['Species'] = mimag_df['Classification'].str.extract(r's__([^;]+)')
base_df = mimag_df[['Bin ID', 'Species']].copy()

# Filter pathogens
pathogens = [
    'Campylobacter_D upsaliensis', 'Enterobacter hormaechei_A', 'Enterobacter roggenkampii',
    'Enterococcus faecalis', 'Enterococcus_B faecium', 'Enterococcus_B hirae',
    'Enterococcus_B lactis', 'Enterococcus_D gallinarum', 'Helicobacter_A bilis_A',
    'Helicobacter_A rappini', 'Helicobacter_B canis', 'Helicobacter_B canis_B',
    'Helicobacter_C magdeburgensis', 'Klebsiella pneumoniae', 'Proteus mirabilis'
]
base_df = base_df[base_df['Species'].isin(pathogens)]

# Count ARGs per Bin
arg_counts = args_df.pivot_table(index='Bin ID', columns='Best_Hit_ARO', aggfunc='size', fill_value=0).reset_index()

# Merge data
full_df = pd.merge(base_df, arg_counts, on='Bin ID', how='left').fillna(0)

# Convert ARG columns to integers
for col in full_df.columns[2:]:
    full_df[col] = full_df[col].astype(int)

# Collapse by species, calculate median
collapsed_df = full_df.groupby('Species').median(numeric_only=True).reset_index()

# Add total ARGs column
collapsed_df['Total_ARGs'] = collapsed_df.iloc[:, 1:].sum(axis=1)
print(collapsed_df)

# We need to finalise the important ARGs for the heatmap-pending
IN2CM = 2.54  
fig, ax = plt.subplots(figsize=(17.2 / IN2CM, 17.2 / IN2CM))

# Filter non-zero ARGs
arg_columns = collapsed_df.columns[1:-1]
non_zero_args = arg_columns[collapsed_df[arg_columns].sum() > 0]
filtered_df = collapsed_df[['Species'] + list(non_zero_args)]
heatmap_df = filtered_df.set_index('Species')
heatmap_df_log = np.log1p(heatmap_df)

# Plot heatmap
sns.heatmap(
    heatmap_df_log,
    cmap='viridis',
    annot=False,
    cbar_kws={'shrink': 0.4, 'aspect': 21, 'pad': 0.01},
    xticklabels=True,
    yticklabels=True,
    square=True,
    ax=ax
)
# Configure colorbar
cbar = ax.collections[0].colorbar
cbar.ax.set_ylabel('')
cbar.ax.tick_params(labelsize=6)

# Configure axes
ax.tick_params(axis='x', rotation=90, labelsize=6)
ax.tick_params(axis='y', labelsize=6)

fig.tight_layout()
fig.savefig('heatmap_args.png', bbox_inches='tight')
fig.show()
