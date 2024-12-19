# -*- coding: utf-8 -*-
"""
Created on Tue 26 Nov 2024 10:02:37
@author: Anna Cusco
"""

import os
import pandas as pd

os.chdir('/data/Projects/ShanghaiDogs')

# Import singleM marker genes
singlem_markers = pd.read_csv('intermediate-outputs/singlem_profiling/beta-div/singlem_metapackage_describe.txt', \
                              sep='\t', encoding= 'unicode_escape', index_col=0)
otu_marker_ls = singlem_markers.index.to_list()
print(otu_marker_ls[0])

# Import metadata
metadata = pd.read_csv('data/ShanghaiDogsMetadata/REP_canid_metadata.csv', \
                       sep=',', encoding= 'unicode_escape', index_col=0)
# n=294 representative samples

# OTU tables filtering functions

def filt_low_OTU(OTUs_tab, min_total, min_prev):
    """
    filt_low_OTU filters out OTUs from an OTU table that do not have
    a min OTU total sum (min_total) or a min prevalence (0-1, min_prev).
    It also prints out removed OTUs. Columns should contain OTU IDs.
    """
    print(OTUs_tab.columns[0]+': this should be a OTU ID')
    otus_count_filt = OTUs_tab.loc[:, OTUs_tab.sum(axis=0) >= min_total]
    prevalence = (otus_count_filt > 0).sum(axis=0) / OTUs_tab.shape[0]
    otus_filt = otus_count_filt.loc[:, prevalence >= min_prev]
    otus_filt_LOW = OTUs_tab.loc[:, ~OTUs_tab.columns.isin(otus_filt.columns)]
    return otus_filt, otus_filt_LOW

def filt_samples_low_OTU(OTUs_tab,min):
    """
    filt_samples_low_OTU filters out samples from an OTU table
    that do not have a min OTU value. It also prints out removed items.
    Index (rows) should contain the samples
    """
    print(OTUs_tab.index[0]+': this should be a sample ID')
    otus_filt = OTUs_tab[OTUs_tab.sum(axis=1) >= min]
    otus_filt_LOW = OTUs_tab[OTUs_tab.sum(axis=1) < min].index.to_frame()
    otus_filt_LOW = pd.merge(otus_filt_LOW,metadata,right_index=True,left_index=True)
    return otus_filt, otus_filt_LOW

def remove_0_sum(OTUs_tab):
    """
    remove_0_sum filters out rows and columns
    from a OTU table that add up to 0. It outputs
    filtered mOTU_tab and a list of removed items.
    """
    otus_filt_wo0 = OTUs_tab[OTUs_tab.sum(axis=1) > 0] #any
    otus_filt_wo0 = otus_filt_wo0.loc[:,otus_filt_wo0.sum(axis=0) > 0]
    otus_filt_0 = OTUs_tab[OTUs_tab.sum(axis=1) == 0].index.tolist()
    otus_filt_0_col = OTUs_tab.loc[:,OTUs_tab.sum(axis=0) == 0]
    otus_filt_0_col = list(otus_filt_0_col)
    removed_0 = otus_filt_0 + otus_filt_0_col
    return otus_filt_wo0, removed_0

# Import OTU tables per each marker gene and check total OTU count per sample in the representative set of samples
summary_otus = pd.DataFrame(index=otu_marker_ls)

for m in otu_marker_ls:
    otu_marker = m
    input_path = ('intermediate-outputs/singlem_profiling/beta-div/unifrac-otu/all-otu-table.'+otu_marker+'.ebd')
    otu_tab = pd.read_csv(input_path, sep='\t', index_col=0)

    ## 1) Reformat OTU table
    # Reformat index names: remove suffixes
    otu_tab.index = otu_tab.index.str.replace('_350', '')
    otu_tab.index = otu_tab.index.str.replace('_1', '')
    otu_tab.index = otu_tab.index.str.replace('_EKDN.*', '', regex=True)
    # Reformat index names: rename run_ids to Biosample_ids
    SRA_metadata = pd.read_csv('external-data/data/dog_microbiome_archive_otu_tables/SRR_Acc_List_Metadata.txt',
                               sep='|', skiprows=[1], index_col=0)
    run_to_biosample = SRA_metadata[[' biosample      ']].reset_index()
    run_to_biosample.columns = ['run', 'biosample']
    run_to_biosample['run'] = run_to_biosample['run'].str.replace(' ', '', regex=True)
    run_to_biosample['biosample'] = run_to_biosample['biosample'].str.replace(' ', '', regex=True)
    run_to_biosample = run_to_biosample.set_index('run')
    for idx in otu_tab.index:
        if (idx.startswith('SRR') or idx.startswith('ERR')) and idx in run_to_biosample.index:
            otu_tab.rename(index={idx: run_to_biosample.loc[idx, 'biosample']}, inplace=True)

    ## 2) Filtering OTU table
    otus_filt_0, otus_filt_rm_otus = filt_low_OTU(otu_tab, 10, 0.02) # filt out OTUs
    otus_filt_1, otus_filt_rm_samples = filt_samples_low_OTU(otus_filt_0, 100) # filt out samples
    otus_filt, rm_all_items = remove_0_sum(otus_filt_1) # filt rows/columns with 0s
    otus_filt.drop('D024', axis=0, inplace=True, errors='ignore')  # filt D024 - in ATBs treatment

    ## 3) Merge with metadata
    otus_w_metadata = pd.merge(otus_filt, metadata[['Study']], left_index=True, right_index=True)

    ## 4) Calculate required statistics
    num_samples_raw = len(otu_tab) # Count of diff samples in otu_tab
    num_otus_raw = len(otu_tab.columns) # Count of diff OTUs in otu_tab
    num_rep_samples = len(otus_w_metadata) # Num representative samples after QC filtering
    num_otus_filt = len(otus_w_metadata.columns)-1 # total filt diff OTU
    filt_otu_count = otus_w_metadata.iloc[:, :-1].sum().sum()  # Sum of all counts in otu_tab
    num_studies = len(otus_w_metadata['Study'].unique())

    # Append information to the DataFrame
    summary_otus.loc[otu_marker, 'initial_num_samples'] = num_samples_raw
    summary_otus.loc[otu_marker, 'initial_num_otus'] = num_otus_raw
    summary_otus.loc[otu_marker, 'filt_num_samples (REP)'] = num_rep_samples
    summary_otus.loc[otu_marker, 'filt_num_otus'] = num_otus_filt
    summary_otus.loc[otu_marker, 'filt_total_sum_otus'] = filt_otu_count
    summary_otus.loc[otu_marker, 'num_studies_final'] = num_studies

# Print and save final table

summary_otus_df = pd.merge(summary_otus,singlem_markers['target_domains'],right_index=True,left_index=True)
summary_otus_df.to_csv('intermediate-outputs/singlem_profiling/beta-div/summary_otus_count.csv',sep=',')
summary_otus_df = pd.read_csv('intermediate-outputs/singlem_profiling/beta-div/summary_otus_count.csv'
                              ,sep=',',index_col=0)

# We decide to use all the markers that target Bacteria and Archaea, and compute median alpha and beta div values
