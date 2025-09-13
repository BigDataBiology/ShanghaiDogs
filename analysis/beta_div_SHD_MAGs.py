# -*- coding: utf-8 -*-
# date "+Created on %a %d %b %Y %H:%M:%S" >> 1-Beta_diversity.py

"""
Created on Thu 06 Jun 2024 16:24:31
@author: Anna Cusco
"""

import pandas as pd
import numpy as np
import skbio
import matplotlib.pyplot as plt

## Import input files

# Import singleM abundance table
singlem_RA = pd.read_csv('intermediate-outputs/singlem_profiling/otus_tab/ALL_OTU_RA_filt_S3.5.rib_prot_S2_rpsB.csv',index_col=0)
singlem_RA = singlem_RA[singlem_RA.index.str.contains('D0')]

# Import MAG abundance table
MAG_RA = pd.read_csv('intermediate-outputs/repbin_coverage_rmean.tsv', sep='\t',index_col=0)
text_to_remove = '_SR_to_95_ANI'
MAG_RA.columns = [col.replace(text_to_remove, '') for col in MAG_RA.columns]

# Import SHD metadata table
metadata = pd.read_csv('data/ShanghaiDogsMetadata/SH_Dog_metadata_red.csv', sep=',', index_col=0)
metadata_T = metadata.T
metadata_T['D000']=metadata_T['D008']
metadata = metadata_T.T

# Import MIMAG table
MIMAG_tab = pd.read_csv('data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv', sep=',')
MIMAG_tab['Bin ID'] = MIMAG_tab['Bin ID'].str.split('.').str[0]

# Remove features with low mean relative abundance
MAG_RA = MAG_RA
singlem_RA = singlem_RA.T

def rm_low_abd_features(abd_tab,min):
    """
    Feature ids need to be on the rows!
    calculate mean mOTU relative abundance and remove those
    motus with a lower mean abundance (min).
    """
    print(abd_tab.index + ': these should be the feature IDs')
    abd_tab['Mean'] = abd_tab.mean(axis=1)
    abd_tab_filt = abd_tab.loc[abd_tab['Mean'].between(min, 1)]
    abd_tab_filt = abd_tab_filt.drop(['Mean'], axis=1)
    return abd_tab_filt

MAG_tab_filt = rm_low_abd_features(MAG_RA,0.0001)
singlem_tab_filt = rm_low_abd_features(singlem_RA,0.0001)

## Log transformation of the RA_abd_tables
def abd_tab_transform(abd_tab,scale_factor):
    """
    Multiply relative abundances by scale_factor. Calculate a
    pseudocount by taking the minimum value from mOTU table
    and divide it by 2. Substitute zeros for the pseudocount.
    Perform log10 transformation of the matrix.

    Required in abd_tab: feature IDs need to be on the columns.
    """
    abd_tab = abd_tab * scale_factor
    np.set_printoptions(precision=5)
    abd_tab_np = abd_tab.to_numpy()

    minval = np.min(abd_tab_np[np.nonzero(abd_tab_np)])
    pseudo = minval/2
    abd_tab.replace([0],pseudo, inplace=True)

    abd_tab = abd_tab.T
    np.set_printoptions(precision=5)
    tab_rows = abd_tab.index.tolist()
    tab_cols = abd_tab.columns.tolist()
    abd_tab_np = abd_tab.to_numpy()

    log_abd_tab = np.log10(abd_tab_np)
    abd_tab_df = pd.DataFrame(log_abd_tab, columns=tab_cols, index=tab_rows)

    return log_abd_tab, tab_rows, tab_cols, abd_tab_df

log_MAG_tab, MAG_rows, MAG_cols, MAG_df = abd_tab_transform(MAG_tab_filt,100000000000)
log_singlem_tab, singlem_rows, singlem_cols, singlem_df = abd_tab_transform(singlem_tab_filt,100000000000)


### COMPUTE AND PLOT BETA DIVERSITY - SHANGHAI PET DOGS
variable = 'Age_classification'
title = 'Sample colored by ' + variable
table = log_MAG_tab
ids = MAG_rows

# Simplify metadata for plotting
mask_unique_var = metadata[variable].value_counts() == 1
unique_var = mask_unique_var[mask_unique_var].index
metadata[variable] = metadata[variable].apply(lambda x: 'n=1' if x in unique_var else x)

bc_div = skbio.diversity.beta_diversity('braycurtis', table, ids=ids, validate=True)
bc_pcoa = skbio.stats.ordination.pcoa(bc_div)

fig, ax = plt.subplots(figsize=(12, 12))
bc_pcoa.plot(metadata,variable,axis_labels=('PC 1', 'PC 2', 'PC 3'),
             title=title, cmap='tab10', s=24)
plt.tight_layout()
plt.show()
