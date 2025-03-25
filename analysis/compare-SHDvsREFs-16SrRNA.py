# -*- coding: utf-8 -*-

"""
Created on Fri Apr 12 12:07:54 2024
@author: Anna Cusco
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import shapiro

os.chdir('/data/Projects/ShanghaiDogs/')
plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

# Import quality reports
SHD_qual = pd.read_csv('data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv',sep=',')
GTDB_qual = pd.read_csv('external-data/data/NCBI_genomes_ref/NCBI_genomes_qual_MIMAG_report.csv',sep=',')

# Filter Representative genomes from Shanghai dogs
SHD_qual_rep = SHD_qual.query('Representative=="Yes" and Quality=="high-quality"')

# Merge the the two qual_report
merged = pd.merge(SHD_qual_rep,GTDB_qual,left_on='GTDBtk fastani Ref',right_on='Name') # include only those that are shared
merged_hq_only = merged[merged['Quality_y'].str.contains('high-quality')]

### Scatter plots for 16S rRNA: RefSeq/GenBank vs HQ MAGs
merged_hq_only = merged_hq_only.sort_values(by='16S rRNA_x', ascending=True)
merged_hq_only['Reference'] = merged_hq_only['GTDBtk fastani Ref'].str.startswith('GCA_').map({True: 'GenBank', False: 'RefSeq'})

# Filter data for GenBank and RefSeq
genbank_data = merged_hq_only[merged_hq_only['Reference'] == 'GenBank']
refseq_data = merged_hq_only[merged_hq_only['Reference'] == 'RefSeq']

## GenBank scatter plot
width_mm = 100
height_mm = 45
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)
ax.scatter(genbank_data['Classification'], genbank_data['16S rRNA_y'],
           label='GenBank genomes', alpha=1, s=6, c='#e6ab02')
ax.scatter(genbank_data['Classification'], genbank_data['16S rRNA_x'],
           label='Shanghai Dogs', alpha=0.3, s=20, c='#1b9e77')
ax.set_xlabel('Species')
ax.set_title('Number of 16S rRNAs',size=10)
ax.set_xticks([])
ax.set_ylim(-0.3, 16)
ax.set_yticks(np.arange(0, 17, 4))
ax.legend()
#sns.despine()
plt.tight_layout()
#plt.show()
fig.savefig('intermediate-outputs/figures/GenBank-vs-SHD_16S_scatterplot.svg')

## RefSeq scatter plot
width_mm = 140
height_mm = 55
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)
ax.scatter(refseq_data['Classification'], refseq_data['16S rRNA_y'],
           label='RefSeq genomes', alpha=1, s=6, c='#a6761d')
ax.scatter(refseq_data['Classification'], refseq_data['16S rRNA_x'],
           label='Shanghai Dogs', alpha=0.3, s=20, c='#1b9e77')
ax.set_xlabel('Species')
ax.set_title('Number of 16S rRNAs',size=10)
ax.set_xticks([])
ax.set_ylim(-0.3, 19)
ax.set_yticks(np.arange(0, 19, 4))
ax.legend()
sns.despine()
plt.tight_layout()
#plt.show()
fig.savefig('intermediate-outputs/figures/RefSeq-vs-SHD_16S_scatterplot.svg')

### Check statistical differences in the number of ribosomal genes
## 1) Assess normality - reshape to long format + Shapiro test
compare_df = refseq_data #genbank_data refseq_data
df = compare_df[['Bin ID','16S rRNA_x','16S rRNA_y']]
df = df.melt(id_vars="Bin ID", var_name="Condition", value_name='#16S rRNA')

conditions = df["Condition"].unique()
normality_results = {}
for condition in conditions:
    stat, p = shapiro(df[df["Condition"] == condition]['#16S rRNA'])
    normality_results[condition] = p
    print(f"Shapiro-Wilk test for {condition}: p-value = {p:.4f}")

## 2) Compute Wilcoxon test for SHD vs REFs 16S rRNA
# Done at compare-SHDvsREFs-contiguity-tRNAs.py