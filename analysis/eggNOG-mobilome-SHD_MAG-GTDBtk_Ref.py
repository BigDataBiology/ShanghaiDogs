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

# Import eggNOG_annotation
# count COG categories

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

# Count 'hits' within each technique by COG category (use only Rep genomes from SHD)
Rep_COGs = pd.merge(all_COGs_descript,SHD_qual[['Bin ID','Representative']],left_on='Name_x',right_on='Bin ID')
Rep_COGs = Rep_COGs.query('Representative == "Yes"')

mobilome_counts_COG = Rep_COGs.groupby('COG_id')[['Count_x','Count_y']].sum().reset_index()
mobilome_counts_COG.columns = ['COG_id','SHD Counts','GTDBtk Ref Counts']
mobilome_counts_COG['prop']=mobilome_counts_COG['GTDBtk Ref Counts']/mobilome_counts_COG['SHD Counts']
mobilome_counts_COG = pd.merge(mobilome_counts_COG,COG_X[['COG','Annotation']],left_on='COG_id',right_on='COG')
mobilome_counts_COG.drop('COG',axis=1,inplace=True)

# Boxplot
sub_mob_counts_COG = mobilome_counts_COG[mobilome_counts_COG['SHD Counts']>=200]

fig, ax = plt.subplots()
ax.clear()
sns.boxplot(data=[sub_mob_counts_COG['SHD Counts'], sub_mob_counts_COG['GTDBtk Ref Counts']], palette="Dark2",ax=ax)
plt.ylabel('Mobilome COGs (hits)')
plt.xticks(ticks=[0, 1], labels=['SHD MAGs (here)','Ref genomes'])
sns.despine(fig, trim=False)
plt.show()

#fig.savefig('analysis/figures/boxplot_COGs_mobilome.svg')

# Heatmap
sub_mob_counts_COG = sub_mob_counts_COG.set_index('COG_id')

fig,ax = plt.subplots(figsize=(7,9))
sns.heatmap(
        sub_mob_counts_COG[['SHD Counts', 'GTDBtk Ref Counts']],
            cmap='Blues', annot=True, fmt="g", linewidths=0.5,linecolor='black')
plt.xlabel('Count')
plt.ylabel('Name')
plt.tight_layout()
plt.show()

#fig.savefig('analysis/figures/heatmap_COGs.svg')

# Add taxonomic information to all_COGs and assess mobilome info by species
all_COGs_w_tax = pd.merge(all_COGs_descript,SHD_qual[['Bin ID','Classification','MIMAG']],left_on='Name_x',right_on='Bin ID')
all_COGs_MIMAG = all_COGs_w_tax.query('MIMAG == "Yes"')
all_COGs_MIMAG.drop(['Bin ID','MIMAG'],axis=1,inplace=True)
all_COGs_MIMAG = all_COGs_MIMAG[all_COGs_MIMAG['GTDBtk fastani Ref'].isin(GTDB_MIMAG_ls)]

# Performed with genomes that fulfill MIMAG criteria (ref, and SHD)
SHD_MIMAG = SHD_qual.query('MIMAG == "Yes"')
count_sp = SHD_MIMAG[('Classification')].value_counts().reset_index()
abd_sp = count_sp.query('count > 30')
abd_sp_ls = list(abd_sp['Classification'])

# Get the mobilome counts for SHD genomes
all_COGs_MIMAG_filt = all_COGs_MIMAG[all_COGs_MIMAG['Classification'].isin(abd_sp_ls)]
SHD_COGs_count = all_COGs_MIMAG_filt[['Name_x','Count_x']]
SHD_COGs_count = SHD_COGs_count.groupby('Name_x').sum()
SHD_COGs_count = pd.merge(SHD_COGs_count,SHD_qual[['Bin ID','Classification']],left_index=True,right_on='Bin ID',how='left')
SHD_COGs_count['Species'] = SHD_COGs_count['Classification'].str.split(';').str[-1].str.split('__').str[-1]
SHD_COGs_ls = list(SHD_COGs_count['Classification'])

# Get the mobilome counts for Ref_GTDB genomes
ext_COGs_count = ext_mobilome_df[['Name','Count']]
ext_COGs_count = ext_COGs_count.groupby('Name').sum()
ext_COGs_count = pd.merge(ext_COGs_count,SHD_qual[['GTDBtk fastani Ref','Classification']],
                          left_index=True,right_on='GTDBtk fastani Ref',how='left')
ext_COGs_count = ext_COGs_count.drop_duplicates()
ext_COGs_count['Species'] = ext_COGs_count['Classification'].str.split(';').str[-1].str.split('__').str[-1]
ext_COGs_count_filt = ext_COGs_count[ext_COGs_count['Classification'].isin(SHD_COGs_ls)]

# check contig number of Ref_GTDB genomes
check = pd.merge(ext_COGs_count_filt,GTDB_qual[['Name','Number']],left_on='GTDBtk fastani Ref',right_on='Name')

# Plot Boxplot + add ext_ref dots

fig,ax = plt.subplots(figsize=(5,10))
sns.boxplot(x=SHD_COGs_count['Species'],y=SHD_COGs_count['Count_x'],ax=ax)
sns.stripplot(x=SHD_COGs_count['Species'], y=SHD_COGs_count['Count_x'], color='black', size=4, ax=ax)
sns.stripplot(x=ext_COGs_count_filt['Species'], y=ext_COGs_count_filt['Count'], color='red', size=4, ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=30,ha='right')
ax.set_xlabel('Species (>30 SHD MAGs)')
ax.set_ylabel('Hits to Mobilome (COG IDs)')
plt.tight_layout()
plt.show()

### Check 'Mobilome' of the MIMAG Ref Genomes that have >50 contigs

GTDB_MIMAG_multiple_contigs = GTDB_MIMAG.query('Number > 50')
GTDB_MIMAG_multiple_contigs_ls = list(GTDB_MIMAG_multiple_contigs['Name'])
abd_sp = count_sp.query('count > 10')
abd_sp_ls = list(abd_sp['Classification'])

# Get the mobilome counts for SHD genomes
all_COGs_MIMAG_filt = all_COGs_MIMAG[all_COGs_MIMAG['GTDBtk fastani Ref'].isin(GTDB_MIMAG_multiple_contigs_ls)]
all_COGs_MIMAG_filt = all_COGs_MIMAG_filt[all_COGs_MIMAG['Classification'].isin(abd_sp_ls)]
SHD_COGs_count = all_COGs_MIMAG_filt[['Name_x','Count_x']]
SHD_COGs_count = SHD_COGs_count.groupby('Name_x').sum()
SHD_COGs_count = pd.merge(SHD_COGs_count,SHD_qual[['Bin ID','Classification']],left_index=True,right_on='Bin ID',how='left')
SHD_COGs_count['Species'] = SHD_COGs_count['Classification'].str.split(';').str[-1].str.split('__').str[-1]
SHD_COGs_ls = list(set(SHD_COGs_count['Classification']))

# Get the mobilome counts for Ref_GTDB genomes
ext_COGs_count_filt = ext_COGs_count[ext_COGs_count['Classification'].isin(SHD_COGs_ls)]

# Plot Boxplot + add ext_ref dots

fig,ax = plt.subplots(figsize=(5,10))
sns.boxplot(x=SHD_COGs_count['Species'],y=SHD_COGs_count['Count_x'],ax=ax)
sns.stripplot(x=SHD_COGs_count['Species'], y=SHD_COGs_count['Count_x'], color='black', size=4, ax=ax)
sns.stripplot(x=ext_COGs_count_filt['Species'], y=ext_COGs_count_filt['Count'], color='red', size=4, ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=30,ha='right')
ax.set_xlabel('Species (>5 SHD MAGs)')
ax.set_ylabel('Hits to Mobilome (COG IDs)')
plt.tight_layout()
plt.show()
