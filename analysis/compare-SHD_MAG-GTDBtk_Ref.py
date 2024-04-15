# -*- coding: utf-8 -*-

"""
Created on Fri Apr 12 12:07:54 2024
@author: Anna Cusco
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

os.chdir('/data/Projects/ShanghaiDogs/')

# Import quality reports
SHD_qual = pd.read_csv('data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv',sep=',')
GTDB_qual = pd.read_csv('external-data/data/NCBI_genomes_ref/NCBI_genomes_qual_MIMAG_report.csv',sep=',')

# Filter Representative genomes from Shanghai dogs
SHD_qual_rep = SHD_qual[SHD_qual['Representative']=='Yes']
SHD_qual_rep = SHD_qual_rep[SHD_qual_rep['Quality']=='high-quality']

# Merge the the two qual_report
merged = pd.merge(SHD_qual_rep,GTDB_qual,left_on='GTDBtk fastani Ref',right_on='Name') # include only those that are shared

# Plots
# Scatterplot
y_SHD = merged['Nr contigs']
x_SHD = merged['16S total_x']
y_Ref = merged['Number']
x_Ref = merged['16S total_y']

plt.scatter(x_Ref, y_Ref, label='Genbank/RefSeq genomes',alpha=0.5)
plt.scatter(x_SHD, y_SHD, label='Shanghai Dogs',alpha=0.5)
plt.xlabel('16S genes')
plt.ylabel('Nr contigs')
plt.title('Scatterplot of Contiguity vs 16S genes')
plt.legend()
#plt.show()
plt.savefig('analysis/figures/SHDvsRef_16S_scatterplot.svg')

# Boxplot
fig, ax = plt.subplots()
ax.clear()
sns.boxplot(data=[x_Ref, x_SHD], palette="Dark2",ax=ax)
plt.ylabel('Nr 16S genes')
plt.xticks(ticks=[0, 1], labels=['Ref genomes', 'SHD MAGs (here)'])
sns.despine(fig, trim=False)
#plt.show()

fig.savefig('analysis/figures/SHDvsRef_16S_boxplot.svg')

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
fig.savefig('analysis/figures/SHDvsRef_16S_histogram.svg')

# Compare the genomes/MAGs that pass MIMAG criteria
GTDB_qual_MIMAG = GTDB_qual[GTDB_qual['MIMAG']=='Yes']
GTDB_qual_MIMAG = GTDB_qual_MIMAG[GTDB_qual_MIMAG['Name'].str.contains('GCF')] # only Ref_seq genome assemblies
SHD_qual_MIMAG = SHD_qual_rep[SHD_qual_rep['MIMAG']=='Yes']
merged_MIMAG = pd.merge(SHD_qual_MIMAG,GTDB_qual_MIMAG,left_on='GTDBtk fastani Ref',right_on='Name') # only shared

# Evaluate the number of ribosomal genes
# First col will store SHD values, second col RefSeq values

counts = {'same_count_ribosomals': [0, 0],
          'two_count_ribosomals': [0, 0],
          'diff_count_ribosomals': [0, 0]}

# Iterate over each row of the dataframe
for index, row in merged_MIMAG.iterrows():
    if row['16S total_x'] == row['23S total_x'] == row['5S total_x']:
        counts['same_count_ribosomals'][0] += 1
    elif row['16S total_x'] == row['23S total_x'] or row['16S total_x'] == row['5S total_x'] or row['23S total_x'] == row['23S total_x']:
        counts['two_count_ribosomals'][0] += 1
    else:
        counts['diff_count_ribosomals'][0] += 1
    if row['16S total_y'] == row['23S total_y'] == row['5S total_y']:
        counts['same_count_ribosomals'][1] += 1
    elif row['16S total_y'] == row['23S total_y'] or row['16S total_y'] == row['5S total_y'] or row['23S total_y'] == row['23S total_y']:
        counts['two_count_ribosomals'][1] += 1
    else:
        counts['diff_count_ribosomals'][1] += 1

ribosomals = pd.DataFrame(counts, index=['ShanghaiDogs: 16S-23S-5S', 'RefSeq: 16S-23S-5S'])

# Compare 16S rRNA genes
df_16S =  merged_MIMAG[['Classification','16S total_x','16S total_y','23S total_x','23S total_y','Nr contigs','Number']]
df_16S['Dif_16S']=df_16S['16S total_x'] - df_16S['16S total_y']
df_16S['Dif_23S']=df_16S['23S total_x'] - df_16S['23S total_y']
df_16S['Dif_Sum']=df_16S['Dif_16S']+df_16S['Dif_23S']
df_16S = df_16S.sort_values(by='Dif_16S')
df_16S = df_16S.drop(columns=['Dif_Sum'])

df_16S['Species']=df_16S['Classification'].str.replace(r'^.*s__', '', regex=True)

# Plotting
labels = df_16S['Species']
values = df_16S['Dif_16S']
plt.figure(figsize=(8, 20))
bars = plt.barh(labels, values, color=['green' if value > 0 else 'red' for value in values],edgecolor='black',linewidth=0.5)

contiguity_ref = df_16S['Number']
for bar, value in zip(bars, contiguity):
    if value == 1:
        bar_position = bar.get_y() + bar.get_height() / 2
        plt.plot(-6, bar_position, marker='_', markersize=10, color='black', linestyle='None')

contiguity_SHD = df_16S['Nr contigs']
for bar, value in zip(bars, contiguity_SHD):
    if value == 1:
        bar_position = bar.get_y() + bar.get_height() / 2
        plt.plot(-5, bar_position, marker='_', markersize=10, color='black', linestyle='None')

plt.xlabel('Dif in N# of 16S genes')
plt.title('16S rRNA counts')
plt.gca().invert_yaxis()  # Invert y-axis to have the highest value on top
plt.tight_layout()
plt.show()