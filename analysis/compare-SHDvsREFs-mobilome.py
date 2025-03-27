# -*- coding: utf-8 -*-

"""
Created on Thu Apr 18 08:52:22 2024
@author: Anna Cusco
"""

import os
import pandas as pd
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import shapiro,wilcoxon
from statsmodels.stats.multitest import multipletests

os.chdir('/data/Projects/ShanghaiDogs/')
plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

# Import qual_reports
SHD_qual = pd.read_csv('data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv',sep=',')
SHD_qual['Bin ID'] = SHD_qual['Bin ID'].str.replace('.fna.gz','')
GTDB_qual = pd.read_csv('external-data/data/NCBI_genomes_ref/NCBI_genomes_qual_MIMAG_report.csv',sep=',')

# Make a list of HQ genomes - we only compare high-quality representative genomes between them
SHD_HQ = SHD_qual.query('Quality == "high-quality" and Representative == "Yes"')
SHD_HQ_ls = list(SHD_HQ['Bin ID'])
GTDB_HQ= GTDB_qual.query('Quality == "high-quality"')
GTDB_HQ_ls = list(GTDB_HQ['Name'])

### NO 'X' COG Category in eggNOG: https://github.com/eggnogdb/eggnog-mapper/issues/424
## Import COG_ids for 'COG_category X' according to NCBI
COG_X = pd.read_csv('external-data/data/NCBI_genomes_ref/eggNOG-annot/NCBI_cog_X_table.tsv',sep='\t')
COG_X_ls = list(COG_X['COG'])

## Import eggNOG_annotation and count COG categories
# external REF genomes
ext_path = 'external-data/data/NCBI_genomes_ref/eggNOG-annot/'
ext_mobilome = []

for filename in glob.glob(os.path.join(ext_path, '*.annotations')):
    annot = pd.read_csv(filename,skiprows=range(4),sep='\t')
    pattern = r'(GCA|GCF)_\d+\.\d+'
    match = re.search(pattern, filename)
    name = match.group(0)
    annot['eggNOG_OGs'].fillna('', inplace=True)
    mob_annot = annot[annot['eggNOG_OGs'].str.contains('|'.join(COG_X_ls))]
    mob_counts = mob_annot['eggNOG_OGs'].value_counts().reset_index()
    mob_counts.columns = ['eggNOG_OGs','Count']
    mob_counts['Name'] = name
    ext_mobilome.append(mob_counts)

ext_mobilome_df = pd.concat(ext_mobilome, ignore_index=True)

# SHD MAGs
shd_path = 'intermediate-outputs/eggNOG-annot/'
shd_mobilome = []

for filename in glob.glob(os.path.join(shd_path, '*.annotations')):
    annot = pd.read_csv(filename,skiprows=range(4),sep='\t')
    pattern = r'SHD\d+_\d+'
    match = re.search(pattern, filename)
    name = match.group(0)
    annot['eggNOG_OGs'].fillna('', inplace=True)
    mob_annot = annot[annot['eggNOG_OGs'].str.contains('|'.join(COG_X_ls))]
    mob_counts = mob_annot['eggNOG_OGs'].value_counts().reset_index()
    mob_counts.columns = ['eggNOG_OGs','Count']
    mob_counts['Name'] = name
    shd_mobilome.append(mob_counts)

shd_mobilome_df = pd.concat(shd_mobilome, ignore_index=True)

# Merge REF genomes and SHD MAGs mobilome hits
shd_mobilome_hits = pd.merge(shd_mobilome_df,SHD_qual[['Bin ID','GTDBtk fastani Ref']],left_on='Name',right_on='Bin ID')
shd_mobilome_hits.drop('Name',axis=1,inplace=True)
merged_hits = pd.merge(shd_mobilome_hits,ext_mobilome_df,left_on=['GTDBtk fastani Ref','eggNOG_OGs'],right_on=['Name','eggNOG_OGs'],how='outer')
merged_hits.drop('Name',axis=1,inplace=True)
merged_hits.fillna(0, inplace=True)

# Keep only high-quality genomes (for both MAGs and REFs)
merged_hits = merged_hits[merged_hits['Bin ID'].isin(SHD_HQ_ls)]
merged_hits = merged_hits[merged_hits['GTDBtk fastani Ref'].isin(GTDB_HQ_ls)]

# Extract lower rank COG_ids and check if they are 'mobilome'
def extract_lowest_rank_cogs(entry, mobilome_list):
    parsed_entries = []

    # Parse the input entries
    for e in entry.split(","):
        cog_id, rest = e.split("@")
        rank, taxonomic_level = rest.split("|")
        # Only include entries that contain 'COG'
        if 'COG' in cog_id:
            parsed_entries.append((cog_id, int(rank), taxonomic_level))

    # If no entries with 'COG', return None
    if not parsed_entries:
        return None

    # Find the COG_id with the lowest taxonomic rank (largest number)
    lowest_rank = max([rank for _, rank, _ in parsed_entries])  # Lowest rank has the highest number
    lowest_rank_cogs = [cog_id for cog_id, rank, _ in parsed_entries if rank == lowest_rank]

    # Filter by mobilome_list
    filtered_cogs = [cog_id for cog_id in lowest_rank_cogs if cog_id in mobilome_list]
    return filtered_cogs if filtered_cogs else None

merged_hits['lowest_rank_cogs'] = merged_hits['eggNOG_OGs'].apply(lambda x: extract_lowest_rank_cogs(x, COG_X_ls))

# Merge COG_ID with the annotation
exploded_hits = merged_hits.explode('lowest_rank_cogs', ignore_index=True) #split if multiple hits
exploded_hits = exploded_hits.rename(columns={'lowest_rank_cogs': 'COG_id'})
exploded_hits_detailed = pd.merge(exploded_hits,COG_X[['COG','Annotation']],left_on='COG_id',right_on='COG',how='left')
# Group back into lists if needed
merged_hits_detailed = (
    exploded_hits_detailed
    .groupby(['eggNOG_OGs', 'Count_x', 'Bin ID', 'GTDBtk fastani Ref', 'Count_y'])
    .agg({
        'COG_id': lambda x: list(x.dropna()),       # Keep COG IDs as lists
        'Annotation': lambda x: list(x.dropna())   # Group annotations
    }).reset_index())

# Add taxonomic information to all_COGs and assess mobilome info by species
merged_hits_w_tax = pd.merge(merged_hits_detailed,SHD_qual[['Bin ID','Classification']],left_on='Bin ID',right_on='Bin ID')
merged_hits_w_tax['Species'] = merged_hits_w_tax['Classification'].str.split(';').str[-1].str.replace('s__','')

merged_hits_w_tax = merged_hits_w_tax[['COG_id','Count_x','Count_y','Annotation',
                                             'Species','Bin ID','GTDBtk fastani Ref']]
merged_hits_w_tax.columns = ['COG ID','COG count (MAGs)','COG count (Reference)','Annotation',
                                'Species','Bin ID','Reference genome']

# Save this dataframe as CSV file
merged_hits_w_tax.to_csv("intermediate-outputs/tables/cog-mobilome-hits-HQ.csv", index=False)

### Total Mobilome 'hits' SHD vs REFs
merged_hits_w_tax['Reference']=merged_hits_w_tax['Reference genome']

for n in merged_hits_w_tax.index:
    ref_value = merged_hits_w_tax.loc[n, 'Reference']
    if 'GCF' in ref_value:
        merged_hits_w_tax.loc[n, 'Reference'] = 'RefSeq'
    elif 'GCA' in ref_value:
        merged_hits_w_tax.loc[n, 'Reference'] = 'GenBank'

# Mobilome COG counts by reference genome origin (RefSeq vs GenBank)
mobilome_total = merged_hits_w_tax.groupby(['Bin ID','Reference genome','Reference','Species'])[['COG count (MAGs)', 'COG count (Reference)']].sum()
mobilome_total = mobilome_total.reset_index()

# Save this dataframe as CSV file
mobilome_total.to_csv("intermediate-outputs/tables/species-mobilome-hits-HQ.csv", index=False)

### Total Mobilome 'hits' SHD vs REFs: BOXPLOT
mobilome_melted = pd.melt(mobilome_total, id_vars=['Bin ID','Reference'],
                          value_vars=['COG count (MAGs)', 'COG count (Reference)'],
                          var_name='Genome', value_name='Count')

mobilome_melted['Genome'] = mobilome_melted['Genome'].str.replace('COG count (','')
mobilome_melted['Genome'] = mobilome_melted['Genome'].str.replace(')','')

# Create 'category' column & merge MAGs into a single category for visualization purposes
mobilome_melted['category'] = mobilome_melted['Genome']+'_'+mobilome_melted['Reference']
mobilome_melted['category'] = mobilome_melted['category'].str.replace('MAGs_GenBank','Shanghai Dog MAG')
mobilome_melted['category'] = mobilome_melted['category'].str.replace('MAGs_RefSeq','Shanghai Dog MAG')
mobilome_melted['log count'] = np.log10(mobilome_melted['Count'])

# Plot boxplot
width_mm = 45
height_mm = 46
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

categories = ['Shanghai Dog MAG', 'Reference_RefSeq','Reference_GenBank']
color_palette = {'Shanghai Dog MAG':'#1b9e77','Reference_RefSeq':'#a6761d','Reference_GenBank':'#e6ab02'}

mobilome_melted = mobilome_melted.sort_values(
    by=['category'],
    ascending=[False])

fig, ax = plt.subplots(figsize=figsize_inch)
ax.clear()

sns.boxplot(data=mobilome_melted,
            x='category',y='log count',
            palette=color_palette,
            ax=ax,
            width=0.8,
            linewidth=1,
            fliersize=2)

ax.set_title('')
ax.set_ylabel('log counts')
ax.set_xlabel('')
ax.set_xticklabels([])
ax.tick_params(axis='x', bottom=False)
sns.despine(fig, trim=False)

plt.tight_layout()
plt.ylim(0, 3)
plt.yticks([0, 1, 2])
# plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/intermediate-outputs/figures/sp_MAG-vs-ref_mobilome_boxplot.svg")

### Total Mobilome 'hits' SHD vs REFs: statistical significance
# 1) Assess normality - reshape to long format + Shapiro test
mobilome_melted = pd.melt(mobilome_total, id_vars=['Bin ID','Reference'],
                          value_vars=['COG count (MAGs)', 'COG count (Reference)'],
                          var_name='Genome', value_name='Count')
mobilome_melted['Genome'] = mobilome_melted['Genome'].str.replace('COG count (','')
mobilome_melted['Genome'] = mobilome_melted['Genome'].str.replace(')','')
mobilome_melted['category'] = mobilome_melted['Genome']+'_'+mobilome_melted['Reference']
mobilome_melted['log count'] = np.log10(mobilome_melted['Count'])

categories = mobilome_melted["category"].unique()
print(categories)

normality_results = {}

for category in categories:
    stat, p = shapiro(mobilome_melted[mobilome_melted["category"] == category]['log count'])
    normality_results[category] = p
    print(f"Shapiro-Wilk test for {category}: p-value = {p:.4f}")

## 2) Compute Wilcoxon test for SHD vs GenBank & SHD vs RefSeq
mobilome_total['log count MAG'] = np.log10(mobilome_total['COG count (MAGs)'])
mobilome_total['log count REF'] = np.log10(mobilome_total['COG count (Reference)'])
refseq_data = mobilome_total[mobilome_total['Reference']=='RefSeq']
genbank_data = mobilome_total[mobilome_total['Reference']=='GenBank']

# Perform Wilcoxon test
p_values = []
results = []

stat, p_value = wilcoxon(genbank_data['log count MAG'],genbank_data['log count REF'],alternative='greater')
results.append({'Comparison': 'mobilome MAG vs REFs', 'Statistic': stat, 'p-value': p_value, 'Group': 'GenBank'})
p_values.append(p_value)
stat, p_value = wilcoxon(refseq_data['log count MAG'],refseq_data['log count REF'],alternative='greater')
results.append({'Comparison': 'mobilome MAG vs REFs', 'Statistic': stat, 'p-value': p_value, 'Group': 'RefSeq'})
p_values.append(p_value)

# Correct for multiple testing (e.g., Benjamini-Hochberg)
_, p_values_corrected, _, _ = multipletests(p_values, method='fdr_bh')

# Add corrected p-values to the results DataFrame
for idx, result in enumerate(results):
    result['Corrected p-value'] = p_values_corrected[idx]

# Create DataFrame
comparison_results = pd.DataFrame(results)
comparison_results['Corrected p-value'] = comparison_results['Corrected p-value'].apply(lambda x: f"{x:.2e}")