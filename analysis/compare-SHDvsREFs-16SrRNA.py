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

# Plots
### Scatter plots for RefSeq representatives vs HQ MAGs
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

# Histogram
colors = sns.color_palette("Dark2", 2)
fig, ax = plt.subplots(1, 2, figsize=(10, 5))

# Plot histogram for x_Ref
ax[0].hist(x_Ref, bins=20, alpha=0.5, label='Ref genomes', edgecolor='black',color=colors[0])
ax[0].set_xlabel('Nr 16S genes')
ax[0].set_ylabel('Frequency')
ax[0].set_title('Reference genomes')
ax[0].set_xlim(0,17)

# Plot histogram for x_SHD
ax[1].hist(x_SHD, bins=20, alpha=0.5, label='SHD MAGs (here)', edgecolor='black',color=colors[1])
ax[1].set_xlabel('Nr 16S genes')
ax[1].set_ylabel('Frequency')
ax[1].set_title('Shanghai Dogs MAGs')
ax[1].set_xlim(0,17)

plt.tight_layout()
#plt.show()
fig.savefig('intermediate-outputs/figures/SHDvsRef_16S_histogram.svg')

# Compare the genomes/MAGs that pass MIMAG criteria
GTDB_qual_MIMAG = GTDB_qual.query("MIMAG =='Yes'")
#GTDB_qual_MIMAG = GTDB_qual_MIMAG[GTDB_qual_MIMAG['Name'].str.contains('GCF')] # only Ref_seq genome assemblies
SHD_qual_MIMAG = SHD_qual_rep.query("MIMAG =='Yes'")
merged_MIMAG = pd.merge(SHD_qual_MIMAG,GTDB_qual_MIMAG,left_on='GTDBtk fastani Ref',right_on='Name') # only shared

# Evaluate the number of ribosomal genes
# First col will store SHD values, second col Reference values

counts = {'same_count_ribosomals': [0, 0],
          'two_count_ribosomals': [0, 0],
          'diff_count_ribosomals': [0, 0]}

# Iterate over each row of the dataframe
for index, row in merged_MIMAG.iterrows():
    if row['16S rRNA_x'] == row['23S rRNA_x'] == row['5S rRNA_x']:
        counts['same_count_ribosomals'][0] += 1
    elif row['16S rRNA_x'] == row['23S rRNA_x'] or row['16S rRNA_x'] == row['5S rRNA_x'] or row['23S rRNA_x'] == row['23S rRNA_x']:
        counts['two_count_ribosomals'][0] += 1
    else:
        counts['diff_count_ribosomals'][0] += 1
    if row['16S rRNA_y'] == row['23S rRNA_y'] == row['5S rRNA_y']:
        counts['same_count_ribosomals'][1] += 1
    elif row['16S rRNA_y'] == row['23S rRNA_y'] or row['16S rRNA_y'] == row['5S rRNA_y'] or row['23S rRNA_y'] == row['23S rRNA_y']:
        counts['two_count_ribosomals'][1] += 1
    else:
        counts['diff_count_ribosomals'][1] += 1

ribosomals = pd.DataFrame(counts, index=['ShanghaiDogs: 16S-23S-5S', 'RefSeq: 16S-23S-5S'])

# Compare 16S rRNA genes
df_16S =  merged_MIMAG[['Classification','16S rRNA_x','16S rRNA_y','23S rRNA_x','23S rRNA_y','Nr contigs','Number']]
df_16S['Dif_16S']=df_16S['16S rRNA_x'] - df_16S['16S rRNA_y']
df_16S['Dif_23S']=df_16S['23S rRNA_x'] - df_16S['23S rRNA_y']
df_16S['Dif_Sum']=df_16S['Dif_16S']+df_16S['Dif_23S']
df_16S = df_16S.sort_values(by='Dif_16S')
df_16S = df_16S.drop(columns=['Dif_Sum'])

df_16S['Species']=df_16S['Classification'].str.replace(r'^.*s__', '', regex=True)

# Plotting
labels = df_16S['Species']
values = df_16S['Dif_16S']

fig, ax = plt.subplots(figsize=(8, 20))
bars = ax.barh(labels, values, color=['green' if value > 0 else 'red' for value in values], edgecolor='black', linewidth=0.5)

contiguity_ref = df_16S['Number']
for bar, value in zip(bars, contiguity_ref):
    if value == 1:
        bar_position = bar.get_y() + bar.get_height() / 2
        ax.plot(-6, bar_position, marker='_', markersize=10, color='black', linestyle='None')

contiguity_SHD = df_16S['Nr contigs']
for bar, value in zip(bars, contiguity_SHD):
    if value == 1:
        bar_position = bar.get_y() + bar.get_height() / 2
        ax.plot(-5, bar_position, marker='_', markersize=10, color='black', linestyle='None')

ax.set_xlabel('Dif in N# of 16S genes',size=14)
ax.set_title('16S rRNA counts',size=18)
ax.set_xlim(-7,15)
ax.set_xticks(range(-4, 16, 2))
ax.invert_yaxis()  # Invert y-axis to have the highest value on top
plt.tight_layout()
sns.despine(trim=False)
fig.savefig('intermediate-outputs/figures/16S_diff_count_species.svg')
#plt.show()