# -*- coding: utf-8 -*-
# date "+Created on %a %d %b %Y %H:%M:%S" >> 1-Beta_diversity.py

"""
Created on Sat 09 Jun 2024 10:44:06
@author: Anna Cusco
"""

import os
import pandas as pd
import numpy as np

os.chdir('/data/Projects/ShanghaiDogs/')

# Import otu/MAG table using short-read
MAG_tax = pd.read_csv('intermediate-outputs/MAGs_tax_profiling/SHDs-MAG_RA_tab.csv',index_col=0)
MAG_tax = MAG_tax*100

singleM_tax = pd.read_csv('intermediate-outputs/singlem_profiling/tax-profiles/SHD-tax-profile-species.tsv', \
                      sep='\t', index_col=0)
singleM_tax.columns = singleM_tax.columns.str.replace('_350', '')
singleM_tax.index = singleM_tax.index.to_series().apply(lambda x: x.split(';')[-1])
singleM_tax.index = singleM_tax.index.str.replace(' s__','')

# Import SHD metadata table
metadata = pd.read_csv('data/ShanghaiDogsMetadata/SH_Dog_metadata_red.csv', sep=',', index_col=0)
metadata_T = metadata.T
metadata_T['D000']=metadata_T['D008']
metadata = metadata_T.T

# Compare TOP 10 for each sample
top_MAG = {}
top_SingleM = {}

for col in MAG_tax.columns:
    top_MAG[col] = MAG_tax[col].nlargest(20).index.tolist()

for col in singleM_tax.columns:
    top_SingleM[col] = singleM_tax[col].nlargest(20).index.tolist()

# Compare TOP 10 for each sample
differences = {}

# Iterate over each sample
for col in top_MAG.keys():
    mag_indices = set(top_MAG[col])
    singlem_indices = set(top_SingleM[col])
    mag_only = mag_indices - singlem_indices
    singlem_only = singlem_indices - mag_indices
    common_indices = mag_indices.intersection(singlem_indices)
    differences[col] = {'MAG_only': list(mag_only), 'SingleM_only': list(singlem_only),
                        'Common_indices': list(common_indices)}

# Convert the differences dictionary to a DataFrame
differences_df = pd.DataFrame(differences)
differences_count_df = differences_df.applymap(lambda x: len(x))

# Print it
print(differences_count_df)