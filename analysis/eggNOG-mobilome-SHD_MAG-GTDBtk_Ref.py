# -*- coding: utf-8 -*-

"""
Created on Thu Apr 18 08:52:22 2024
@author: Anna Cusco
"""

import os
import glob
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

os.chdir('/data/Projects/ShanghaiDogs/')
plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

# Import qual_reports
SHD_qual = pd.read_csv('data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv',sep=',')
SHD_qual['Bin ID'] = SHD_qual['Bin ID'].str.replace('.fna.gz','')
SHD_Ref = SHD_qual[['Bin ID','GTDBtk fastani Ref']]
GTDB_qual = pd.read_csv('external-data/data/NCBI_genomes_ref/NCBI_genomes_qual_MIMAG_report.csv',sep=',')
GTDB_MIMAG= GTDB_qual.query('MIMAG == "Yes"')
GTDB_MIMAG_ls = list(GTDB_MIMAG['Name'])

### NO 'X' COG Category in eggNOG!!!!
### https://github.com/eggnogdb/eggnog-mapper/issues/424
# Import COG_ids for COG_category X according to NCBI
COG_X = pd.read_csv('external-data/data/NCBI_genomes_ref/eggNOG-annot/NCBI_cog_X_table.tsv',sep='\t')
COG_X_ls = list(COG_X['COG'])

# Import eggNOG_annotation & count COG categories
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

a = pd.merge(shd_mobilome_df,SHD_Ref,left_on='Name',right_on='Bin ID')
a.drop('Bin ID',axis=1,inplace=True)

all_COGs=pd.merge(a,ext_mobilome_df,left_on=['GTDBtk fastani Ref','eggNOG_OGs'],right_on=['Name','eggNOG_OGs'],how='outer')
all_COGs.fillna(0, inplace=True)

all_COGs['COG_id']=all_COGs['eggNOG_OGs'].str.split('@').str[0]
all_COGs_descript = pd.merge(all_COGs,COG_X[['COG','Annotation']],left_on='COG_id',right_on='COG',how='left')
all_COGs_descript.drop(['COG','eggNOG_OGs','Name_y'],axis=1,inplace=True)
all_COGs_descript = all_COGs_descript.query('`GTDBtk fastani Ref` != "0"')
all_COGs_descript = all_COGs_descript[['COG_id','Name_x','GTDBtk fastani Ref','Count_x','Count_y','Annotation']]

# Count 'hits' within each technique by COG category (use only Representative genomes from SHD)
Rep_COGs = pd.merge(all_COGs_descript,SHD_qual[['Bin ID','Representative']],left_on='Name_x',right_on='Bin ID')
Rep_COGs = Rep_COGs.query('Representative == "Yes"')

mobilome_counts_COG = Rep_COGs.groupby('COG_id')[['Count_x','Count_y']].sum().reset_index()
mobilome_counts_COG.columns = ['COG_id','SHD Counts','GTDBtk Ref Counts']
mobilome_counts_COG['prop']=mobilome_counts_COG['SHD Counts']/mobilome_counts_COG['GTDBtk Ref Counts']
mobilome_counts_COG = pd.merge(mobilome_counts_COG,COG_X[['COG','Annotation']],left_on='COG_id',right_on='COG')
mobilome_counts_COG.drop('COG',axis=1,inplace=True)

# Boxplot
sub_mob_counts_COG = mobilome_counts_COG.query("`SHD Counts` >= 50 and `GTDBtk Ref Counts` >= 25")
sub_mob_counts_COG['prop'].mean()

fig, ax = plt.subplots()
ax.clear()
sns.boxplot(data=[sub_mob_counts_COG['SHD Counts'], sub_mob_counts_COG['GTDBtk Ref Counts']], palette="Dark2",ax=ax)
plt.ylabel('Mobilome COGs (hits)')
plt.xticks(ticks=[0, 1], labels=['SHD MAGs (here)','Ref genomes'])
sns.despine(fig, trim=False)
#plt.show()

fig.savefig('intermediate-outputs/figures/boxplot_total_COGs_mobilome.svg')

# Heatmap
sub_mob_counts_COG['Tag']=sub_mob_counts_COG['COG_id']+'_'+sub_mob_counts_COG['Annotation']
sub_mob_counts_COG = sub_mob_counts_COG.set_index('Tag')
sub_mob_counts_COG = sub_mob_counts_COG.sort_values(by='SHD Counts',ascending=False)

fig,ax = plt.subplots(figsize=(11,9))
sns.heatmap(
    sub_mob_counts_COG[['SHD Counts', 'GTDBtk Ref Counts']],
    cmap='Blues',
    annot=True,
    fmt="g",
    linewidths=0.5,
    linecolor='black',
    ax=ax,
    cbar_kws={'shrink': 0.5})
ax.set_aspect('auto')
plt.ylabel('')
plt.tight_layout()
#plt.show()

fig.savefig('intermediate-outputs/figures/heatmap_specific_COGs_count.svg')

# Add taxonomic information to all_COGs and assess mobilome info by species
all_COGs_w_tax = pd.merge(all_COGs_descript,SHD_qual[['Bin ID','Classification','MIMAG']],left_on='Name_x',right_on='Bin ID')

# Mobilome count for Ref_GTDB genomes
ext_COGs_count = ext_mobilome_df[['Name','Count']]
ext_COGs_count = ext_COGs_count.groupby('Name').sum()
ext_COGs_count = pd.merge(ext_COGs_count,SHD_qual[['GTDBtk fastani Ref','Classification']],
                          left_index=True,right_on='GTDBtk fastani Ref',how='left')
ext_COGs_count = ext_COGs_count.drop_duplicates()
ext_COGs_count['Species'] = ext_COGs_count['Classification'].str.split(';').str[-1].str.split('__').str[-1]

### COMPARE MOBILOME COUNTS SHDs vs REF GENOMES (in high-quality MAGs)
# Filter high-quality and SHD MAGs with 10+ representatives
SHD_HQ = SHD_qual.query('Quality == "high-quality"')
count_sp = SHD_HQ['Classification'].value_counts().reset_index()
abd_sp = count_sp.query('count > 10')
abd_sp_ls = list(abd_sp['Classification'])

# Classify GTDB genomes for contiguity
GTDB_contiguous = GTDB_qual.query('Number < 20')
GTDB_contiguous_ls = list(GTDB_contiguous['Name'])
GTDB_non_contiguous = GTDB_qual.query('Number > 80')
GTDB_non_contiguous_ls = list(GTDB_non_contiguous['Name'])

# Filter all_COGs_w_tax for contiguous Ref genomes and abundant SHD MAGs
all_COGs_w_tax_contiguous = all_COGs_w_tax[all_COGs_w_tax['GTDBtk fastani Ref'].isin(GTDB_contiguous_ls)]
all_COGs_w_tax_contiguous = all_COGs_w_tax_contiguous[all_COGs_w_tax_contiguous['Classification'].isin(abd_sp_ls)]
SHD_COGs_count_contiguous = all_COGs_w_tax_contiguous.groupby('Name_x')['Count_x'].sum().reset_index()
SHD_COGs_count_contiguous = pd.merge(SHD_COGs_count_contiguous, SHD_HQ[['Bin ID', 'Classification']],
                                     left_on='Name_x', right_on='Bin ID', how='left')
SHD_COGs_count_contiguous['Species'] = SHD_COGs_count_contiguous['Classification'].str.split(';').str[-1].str.split('__').str[-1]

# Filter ext_COGs_count for contiguous genomes and SHD abundant species
ext_COGs_count_contiguous = ext_COGs_count[ext_COGs_count['Classification'].isin(SHD_COGs_count_contiguous['Classification'])]
ext_COGs_count_contiguous_MIMAG = ext_COGs_count_contiguous[ext_COGs_count_contiguous['GTDBtk fastani Ref'].isin(GTDB_MIMAG_ls)]

# Filter all_COGs_w_tax for non-contiguous Ref genomes and abundant SHD MAGs
all_COGs_w_tax_non_contiguous = all_COGs_w_tax[all_COGs_w_tax['GTDBtk fastani Ref'].isin(GTDB_non_contiguous_ls)]
all_COGs_w_tax_non_contiguous = all_COGs_w_tax_non_contiguous[all_COGs_w_tax_non_contiguous['Classification'].isin(abd_sp_ls)]
SHD_COGs_count_non_contiguous = all_COGs_w_tax_non_contiguous.groupby('Name_x')['Count_x'].sum().reset_index()
SHD_COGs_count_non_contiguous = pd.merge(SHD_COGs_count_non_contiguous, SHD_HQ[['Bin ID', 'Classification']],
                                         left_on='Name_x', right_on='Bin ID', how='left')
SHD_COGs_count_non_contiguous['Species'] = SHD_COGs_count_non_contiguous['Classification'].str.split(';').str[-1].str.split('__').str[-1]

# Filter ext_COGs_count for non contiguous genomes and SHD abundant species
ext_COGs_count_non_contiguous = ext_COGs_count[ext_COGs_count['Classification'].isin(SHD_COGs_count_non_contiguous['Classification'])]
ext_COGs_count_non_contiguous_MIMAG = ext_COGs_count_non_contiguous[ext_COGs_count_non_contiguous['GTDBtk fastani Ref'].isin(GTDB_MIMAG_ls)]

# Plot Boxplots
fig, axes = plt.subplots(1, 2, figsize=(15, 5))

palette_axes = sns.color_palette("Dark2", 6)
order_cont=SHD_COGs_count_contiguous.groupby('Species')['Count_x'].median().sort_values().index
order_non_cont=SHD_COGs_count_non_contiguous.groupby('Species')['Count_x'].median().sort_values().index

# Plot for contiguous genomes
sns.boxplot(x=SHD_COGs_count_contiguous['Species'], y=SHD_COGs_count_contiguous['Count_x'], ax=axes[0],\
            order=order_cont,color=palette_axes[0],fliersize=0)
sns.stripplot(x=SHD_COGs_count_contiguous['Species'], y=SHD_COGs_count_contiguous['Count_x'], color='black',\
              size=2, ax=axes[0],order=order_cont)
sns.stripplot(x=ext_COGs_count_contiguous['Species'], y=ext_COGs_count_contiguous['Count'], color='red',\
              size=4, ax=axes[0],order=order_cont)
sns.stripplot(x=ext_COGs_count_contiguous_MIMAG['Species'], y=ext_COGs_count_contiguous_MIMAG['Count'], color='blue',\
              size=4, ax=axes[0],order=order_cont)
axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=30, ha='right')
axes[0].set_xlabel('Species (>10 SHD MAGs)')
axes[0].set_ylabel('Hits to Mobilome COGs')
axes[0].set_title('Contiguous Ref Genomes (<20 contigs)')
axes[0].set_ylim(0, 350)

# Plot for many_contigs genomes
sns.boxplot(x=SHD_COGs_count_non_contiguous['Species'], y=SHD_COGs_count_non_contiguous['Count_x'], ax=axes[1],\
            order=order_non_cont,color=palette_axes[1],fliersize=0)
sns.stripplot(x=SHD_COGs_count_non_contiguous['Species'], y=SHD_COGs_count_non_contiguous['Count_x'], color='black',\
              size=2, ax=axes[1],order=order_non_cont)
sns.stripplot(x=ext_COGs_count_non_contiguous['Species'], y=ext_COGs_count_non_contiguous['Count'], color='red',\
              size=4, ax=axes[1],order=order_non_cont)
sns.stripplot(x=ext_COGs_count_non_contiguous_MIMAG['Species'], y=ext_COGs_count_non_contiguous_MIMAG['Count'], color='blue',\
              size=4, ax=axes[1],order=order_non_cont)
axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=30, ha='right')
axes[1].set_xlabel('Species (>10 SHD MAGs)')
axes[1].set_ylabel('Hits to Mobilome COGs')
axes[1].set_title('Non-contiguous Ref Genomes (>80 contigs)')
axes[1].set_ylim(0, 350)

plt.tight_layout()
fig.savefig('intermediate-outputs/figures/Mobilome_by_species.svg')
#plt.show()

### COMPARE MOBILOME COUNTS IN SHDs OF MOST ABUNDANT SPs (ORDERED BY TAXONOMIC CLASSIFICATION)
SHD_MIMAG = SHD_qual.query('MIMAG == "Yes"')
count_sp = SHD_MIMAG['Classification'].value_counts().reset_index()
abd_sp = count_sp.query('count > 20')
abd_sp_ls = list(abd_sp['Classification'])

all_COGs_w_tax_MIMAG = all_COGs_w_tax[all_COGs_w_tax['Classification'].isin(abd_sp_ls)]
all_COGs_w_tax_MIMAG = all_COGs_w_tax_MIMAG.query('MIMAG == "Yes"')
SHD_COGs_count_MIMAG = all_COGs_w_tax_MIMAG.groupby('Name_x')['Count_x'].sum().reset_index()
SHD_COGs_count_MIMAG = pd.merge(SHD_COGs_count_MIMAG, SHD_MIMAG[['Bin ID', 'Classification']],
                                     left_on='Name_x', right_on='Bin ID', how='left')
SHD_COGs_count_MIMAG['Species'] = SHD_COGs_count_MIMAG['Classification'].str.split(';').str[-1].str.replace('s__','')
SHD_COGs_count_MIMAG['Genera'] = SHD_COGs_count_MIMAG['Classification'].str.split(';').str[-2].str.replace('g__','')
SHD_COGs_count_MIMAG['Family'] = SHD_COGs_count_MIMAG['Classification'].str.split(';').str[-3].str.replace('f__','')
SHD_COGs_count_MIMAG['Phylum'] = SHD_COGs_count_MIMAG['Classification'].str.split(';').str[1].str.replace('p__','')

SHD_COGs_count_MIMAG_sorted = SHD_COGs_count_MIMAG.sort_values(by='Count_x')
custom_phylum_order = ['Bacillota','Bacillota_A', 'Bacillota_C', 'Bacteroidota','Pseudomonadota',\
                       'Fusobacteriota','Actinomycetota']
SHD_COGs_count_MIMAG_sorted['Phylum'] = pd.Categorical(SHD_COGs_count_MIMAG_sorted['Phylum'], categories=custom_phylum_order, ordered=True)
SHD_COGs_count_MIMAG_sorted = SHD_COGs_count_MIMAG_sorted.sort_values(by='Phylum')

# Plot boxplot
fig, ax = plt.subplots(figsize=(12, 4))
sns.boxplot(x=SHD_COGs_count_MIMAG_sorted['Species'], y=SHD_COGs_count_MIMAG_sorted['Count_x'], ax=ax,\
            hue=SHD_COGs_count_MIMAG_sorted['Phylum'],palette='Dark2',fliersize=0, width=0.65)

ax.set_xlim(ax.get_xlim()[0] - 0.5, ax.get_xlim()[1] + 0.5)
ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right')
ax.set_xlabel('Species (>20 SHD MAGs)')
ax.set_ylabel('Mobilome hits')
ax.legend(loc='upper left', bbox_to_anchor=(1, 1), title='Phylum')

plt.tight_layout()

fig.savefig('intermediate-outputs/figures/Mobilome_by_species_tax.svg')
#plt.show()

###### Plot Mobilome boxplots by GCA, vs GCF
all_COGs_w_tax_ = pd.merge(all_COGs_w_tax, SHD_qual[['Quality','Bin ID','Representative']],left_on='Bin ID',right_on='Bin ID')
all_COGs_w_tax_ = all_COGs_w_tax_.query('`GTDBtk fastani Ref` != 0 and Representative =="Yes"')
all_COGs_w_tax_['ref_qual']=all_COGs_w_tax_['GTDBtk fastani Ref']

for n in all_COGs_w_tax_.index:
    ref_value = all_COGs_w_tax_.loc[n, 'ref_qual']
    if 'GCF' in ref_value:
        all_COGs_w_tax_.loc[n, 'ref_qual'] = 'RefSeq reference (GCF)'
    elif 'GCA' in ref_value:
        all_COGs_w_tax_.loc[n, 'ref_qual'] = 'GenBank reference (GCA)'

### Mobilome COG counts by reference genome assembly quality

mobilome_MAGs = all_COGs_w_tax_.groupby(['Bin ID','ref_qual'])[['Count_x', 'Count_y']].sum()
mobilome_MAGs = mobilome_MAGs.reset_index()

mobilome_melted = pd.melt(mobilome_MAGs, id_vars=['Bin ID','ref_qual'],
                    value_vars=['Count_x', 'Count_y'],
                    var_name='Genome', value_name='Count')

mobilome_melted['category']=mobilome_melted['Genome']+'_'+mobilome_melted['ref_qual']
mobilome_melted['category']=mobilome_melted['category'].str.replace('Count_x','MAG')
mobilome_melted['category']=mobilome_melted['category'].str.replace('Count_y','REF')

order = ['REF_RefSeq reference (GCF)','MAG_RefSeq reference (GCF)','REF_GenBank reference (GCA)','MAG_GenBank reference (GCA)']
mobilome_melted['category'] = pd.Categorical(mobilome_melted['category'], categories=order, ordered=True)
mobilome_melted = mobilome_melted.sort_values('category')

# If merging canine MAGs (mobilome counts)
mobilome_melted['category'] = mobilome_melted['category'].str.replace('MAG_GenBank reference (GCA)','Shanghai Dog MAG')
mobilome_melted['category'] = mobilome_melted['category'].str.replace('MAG_RefSeq reference (GCF)','Shanghai Dog MAG')
mobilome_melted['log count'] = np.log10(mobilome_melted['Count'])

# Plot boxplot mobilome GCA vs GCF
width_mm = 45
height_mm = 50
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

categories = ['Shanghai Dog MAG','REF_RefSeq reference (GCF)','REF_GenBank reference (GCA)']
color_palette = {'Shanghai Dog MAG':'#1b9e77','REF_RefSeq reference (GCF)':'#a6761d','REF_GenBank reference (GCA)':'#e6ab02'}
mobilome_melted = mobilome_melted.sort_values(
    by=['category'],
    ascending=[False]
)


fig, ax = plt.subplots(figsize=figsize_inch)
ax.clear()
sns.boxplot(data=mobilome_melted,
            x='category',y='log count',
            palette=color_palette,
            ax=ax,
            width=0.7,
            linewidth=1,
            flierprops={
                'marker': 'd',  # Shape of outliers
                'color': 'gray',  # Color of outliers
                'markersize': 2,  # Size of outliers
                'linestyle': 'none'  # No connecting line for outliers
    })
ax.set_ylabel('log10 counts')
ax.set_xlabel('')
ax.set_xticklabels([])
ax.set_title('mobilome')
ax.set_yticks(np.arange(0, 3, 1))
ax.tick_params(axis='x', bottom=False)
sns.despine(fig, trim=False)

plt.tight_layout()
# plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/intermediate-outputs/figures/sp_MAG-vs-ref_mobilome_boxplot.svg")
