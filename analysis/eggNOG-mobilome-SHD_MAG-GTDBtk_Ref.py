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

# Import qual_report
SHD_qual = pd.read_csv('data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv',sep=',')
SHD_qual['Bin ID'] = SHD_qual['Bin ID'].str.replace('.fna.gz','')
SHD_Ref = SHD_qual[['Bin ID','GTDBtk fastani Ref']]


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
  with open(filename, 'r') as file:
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
  with open(filename, 'r') as file:
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
all_COGs.fillna('0', inplace=True)
all_COGs['Count_x']=(all_COGs['Count_x']).astype(float)
all_COGs['Count_y']=(all_COGs['Count_y']).astype(float)

all_COGs['COG_id']=all_COGs['eggNOG_OGs'].str.split('@').str[0]
all_COGs_descript = pd.merge(all_COGs,COG_X[['COG','Annotation']],left_on='COG_id',right_on='COG',how='left')
all_COGs_descript.drop(['COG','eggNOG_OGs','Name_y'],axis=1,inplace=True)
all_COGs_descript = all_COGs_descript[all_COGs_descript['GTDBtk fastani Ref'] != '0']
all_COGs_descript = all_COGs_descript[['COG_id','Name_x','GTDBtk fastani Ref','Count_x','Count_y','Annotation']]

# Count 'hits' within each technique by COG category (use only Rep genomes from SHD)
Rep_COGs = pd.merge(all_COGs_descript,SHD_qual[['Bin ID','Representative']],left_on='Name_x',right_on='Bin ID')
Rep_COGs = Rep_COGs[Rep_COGs['Representative']=='Yes']

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

# Create the heatmap
plt.figure(figsize=(7,9))
sns.heatmap(sub_mob_counts_COG[['SHD Counts', 'GTDBtk Ref Counts']], cmap='Blues', annot=True, fmt="g", linewidths=0.5,linecolor='black')
plt.xlabel('Count')
plt.ylabel('Name')
plt.tight_layout()
plt.show()

#fig.savefig('analysis/figures/heatmap_COGs.svg')
