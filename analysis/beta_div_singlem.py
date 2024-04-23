# -*- coding: utf-8 -*-
"""
Created on Mon 08 Apr 2024 16:16:51
@author: Anna Cusco
"""

import os
import pandas as pd
import numpy as np
import skbio
from skbio.stats.distance import anosim,permanova
import matplotlib.pyplot as plt

os.chdir('/data/Projects/ShanghaiDogs')

# IMPORT INPUT FILES
metadata = pd.read_csv('data/ShanghaiDogsMetadata/ALL_canid_felid_metadata.csv', \
                       sep=',', encoding= 'unicode_escape', index_col=0)
metadata_rep = metadata.query('Genus == "Canis" and Representative=="Yes"')

otu_tab = pd.read_csv('intermediate-outputs/singlem_profiling/beta-div/all-otu-table.S3.5.rib_prot_S2_rpsB.ebd', \
                      sep='\t', index_col=0)
otu_tab.index = otu_tab.index.str.replace('_350','')

# Update the index in the otu_tab to the biosample_id
SRA_metadata = pd.read_csv('external-data/data/dog_microbiome_archive_otu_tables/SRR_Acc_List_Metadata.txt', \
                               sep='|',skiprows=[1],index_col=0)
run_to_biosample = SRA_metadata[[' biosample      ']].reset_index()
run_to_biosample.columns = ['run','biosample']
run_to_biosample['run'] = run_to_biosample['run'].str.replace(' ', '', regex=True)
run_to_biosample['biosample'] = run_to_biosample['biosample'].str.replace(' ', '', regex=True)
run_to_biosample = run_to_biosample.set_index('run')

# rename run_ids to Biosample_ids
for idx in otu_tab.index:
    if idx.startswith('SRR') or idx.startswith('ERR'):
        otu_tab.rename(index={idx: run_to_biosample.loc[idx, 'biosample']}, inplace=True)


# Filter out samples with low OTU count
def filt_samples_low_OTU(OTUs_tab,min):
    """
    filt_samples_low_OTU filters out samples from an OTU table
    that do not have a min OTU value. It also prints out removed items.
    Index (rows) should contain the samples
    """
    print(OTUs_tab.index+': these should be the sample IDs')
    otus_filt = OTUs_tab[OTUs_tab.sum(axis=1) >= min]
    otus_filt_LOW = OTUs_tab[OTUs_tab.sum(axis=1) < min].index.to_frame()
    otus_filt_LOW = pd.merge(otus_filt_LOW,metadata_rep,right_index=True,left_index=True)
    return otus_filt, otus_filt_LOW

otus_all_filt, otus_all_filt_low = filt_samples_low_OTU(otu_tab, 100)

# Remove 0 columns and rows

def remove_0_sum(OTUs_tab):
    """
    remove_0_sum filters out rows and columns
    from a OTU table that add up to 0. It outputs
    filtered mOTU_tab and a list of removed items.
    """
    otus_filt_wo0 = OTUs_tab[OTUs_tab.sum(axis=1) > 0]
    otus_filt_wo0 = otus_filt_wo0.loc[:,otus_filt_wo0.sum(axis=0) > 0]
    otus_filt_0 = OTUs_tab[OTUs_tab.sum(axis=1) == 0].index.tolist()
    otus_filt_0_col = OTUs_tab.loc[:,OTUs_tab.sum(axis=0) == 0]
    otus_filt_0_col = list(otus_filt_0_col)
    removed_0 = otus_filt_0 + otus_filt_0_col
    return otus_filt_wo0, removed_0

otus_all_filt_wo0, rm_all_items = remove_0_sum(otus_all_filt)

# Calculate relative abundance of OTU table
def rel_ab_otu(OTUs_tab):
    """
    check OTU table orientation, tranpose when necessary,
    and calculate relative abundance of OTUs/sample.
    Index (rows) should contain the samples
    """
    OTUs_tab['Sum'] = OTUs_tab.sum(axis=1)
    otus_RA = OTUs_tab.iloc[:, :].div(OTUs_tab.Sum, axis=0)
    otus_RA = otus_RA.drop(['Sum'], axis=1)
    return otus_RA

otus_all_RA = rel_ab_otu(otus_all_filt_wo0)

# Remove OTUs with low mean relative abundance
def rm_low_mean_otus(OTUs_tab,min):
    """
    transpose OTU table (OTU_ids should be on the rows)
    calculate mean mOTU relative abundance and remove those
    motus with a lower mean abundance (min).
    """
    OTUs_tab = OTUs_tab.T
    print(OTUs_tab.index + ': these should be the OTUs')
    OTUs_tab['Mean'] = OTUs_tab.mean(axis=1)
    OTUs_tab_filt = OTUs_tab.loc[OTUs_tab['Mean'].between(min, 1)]
    OTUs_tab_filt = OTUs_tab_filt.drop(['Mean'], axis=1)
    return OTUs_tab_filt

otus_all_RA_filt = rm_low_mean_otus(otus_all_RA,0.00001)

# Filter metadata and otu_tab to contain info for only the samples that are in the PCoA plot
samples_id = list(otus_all_RA_filt.columns)
metadata_filt = metadata_rep[metadata_rep.index.isin(samples_id)]
samples_id = list(metadata_filt.index)
otus_all_RA_filt = otus_all_RA_filt.T
otus_all_RA_filt_2 = otus_all_RA_filt[otus_all_RA_filt.index.isin(samples_id)]

# Log transformation of the OTUs table
def otus_transform(OTUs_tab):
    """
    Multiply relative abundances by 10^6. Calculate a
    pseudocount by taking the minimum value from OTU table
    and divide it by 2. Substitute zeros for the pseudocount.
    Perform log10 transformation of the matrix.

    Required in OTU_tab: OTUs IDs need to be on the columns.
    """
    OTUs_tab = OTUs_tab.T
    print(OTUs_tab.index + ': these should be the OTUs')
    OTUs_tab = OTUs_tab * 1000000
    np.set_printoptions(precision=5)
    otus_np = OTUs_tab.to_numpy()

    minval = np.min(otus_np[np.nonzero(otus_np)])
    pseudo = minval/2
    OTUs_tab.replace([0],pseudo, inplace=True)

    OTUs_tab = OTUs_tab.T
    np.set_printoptions(precision=5)
    otus_rows = OTUs_tab.index.tolist()
    otus_cols = OTUs_tab.columns.tolist()
    otus_np = OTUs_tab.to_numpy()

    log_otus_tab = np.log10(otus_np)
    otus_tab_F = pd.DataFrame(log_otus_tab, columns=otus_cols, index=otus_rows)

    return log_otus_tab, otus_rows, otus_cols, otus_tab_F

log_otus_all_tab, otus_all_rows, otus_all_cols, otus_all_F = otus_transform(otus_all_RA_filt_2)

### COMPUTE AND PLOT BETA DIVERSITY

bc_div = skbio.diversity.beta_diversity('braycurtis', log_otus_all_tab, ids=otus_all_rows, validate=True)
bc_pcoa = skbio.stats.ordination.pcoa(bc_div)
bc_pcoa.plot(metadata_filt,'Study',axis_labels=('PC 1', 'PC 2', 'PC 3'),
             title='Samples colored by Study', cmap='tab20', s=24)
plt.tight_layout()
#plt.show()
plt.savefig('analysis/figures/PCoA_all_study_id.svg')

bc_div = skbio.diversity.beta_diversity('braycurtis', log_otus_all_tab, ids=otus_all_rows, validate=True)
bc_pcoa = skbio.stats.ordination.pcoa(bc_div)
bc_pcoa.plot(metadata_filt,'env_classification',axis_labels=('PC 1', 'PC 2', 'PC 3'),
             title='Samples colored by env_classification', cmap='Dark2', s=24)
plt.tight_layout()
#plt.show()
plt.savefig('analysis/figures/PCoA_all_env_classification.svg')

