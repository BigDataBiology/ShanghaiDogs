# -*- coding: utf-8 -*-
"""
Created on Mon 08 Apr 2024 16:16:51
@author: Anna Cusco
"""

import os
import pandas as pd
import numpy as np
import skbio
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

os.chdir('/data/Projects/ShanghaiDogs')
plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

# IMPORT INPUT FILES

# Metadata
metadata = pd.read_csv('data/ShanghaiDogsMetadata/REP_canid_metadata.csv', \
                       sep=',', encoding= 'unicode_escape', index_col=0)
metadata['env_classification'] = metadata['env_classification'].str.replace(' captive','')

# OTU tables for PCoA plot
otu_marker = 'S3.42.RNA_pol_A_bac' #S3.5.rib_prot_S2_rpsB #S3.1.ribosomal_protein_L2_rplB #S3.42.RNA_pol_A_bac
input_path = ('intermediate-outputs/singlem_profiling/beta-div/unifrac-otu/all-otu-table.'+otu_marker+'.ebd')
otu_tab = pd.read_csv(input_path, sep='\t', index_col=0)

# PRE-PROCESSING OF INPUT FILES

# Reformat index names: remove suffixes
otu_tab.index = otu_tab.index.str.replace('_350','')
otu_tab.index = otu_tab.index.str.replace('_1','')
otu_tab.index = otu_tab.index.str.replace('_EKDN.*','',regex=True)

# Reformat index names: rename run_ids to Biosample_ids
SRA_metadata = pd.read_csv('external-data/data/dog_microbiome_archive_otu_tables/SRR_Acc_List_Metadata.txt', \
                               sep='|',skiprows=[1],index_col=0)
run_to_biosample = SRA_metadata[[' biosample      ']].reset_index()
run_to_biosample.columns = ['run','biosample']
run_to_biosample['run'] = run_to_biosample['run'].str.replace(' ', '', regex=True)
run_to_biosample['biosample'] = run_to_biosample['biosample'].str.replace(' ', '', regex=True)
run_to_biosample = run_to_biosample.set_index('run')

for idx in otu_tab.index:
    if (idx.startswith('SRR') or idx.startswith('ERR')) and idx in run_to_biosample.index:
        otu_tab.rename(index={idx: run_to_biosample.loc[idx, 'biosample']}, inplace=True)

# FILTERING OTU TABLE

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
    otus_filt_LOW = pd.merge(otus_filt_LOW,metadata,right_index=True,left_index=True)
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
#otus_all_RA.to_csv('intermediate-outputs/singlem_profiling/beta-div/unifrac-otu/REP_OTU_RA_filt_S3.5.rib_prot_S2_rpsB.csv')

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

# DEFINE SAMPLES TO INCLUDE ON BETA DIVERSITY
samples_id = list(otus_all_RA_filt.columns) # Removed samples with low mean otus
metadata_filt = metadata[metadata.index.isin(samples_id)]
samples_id = list(metadata_filt.index)

# Further filter out extra samples if necessary
otus_RA = otus_all_RA_filt.T
otus_RA_filt = otus_RA[otus_RA.index.isin(samples_id)]
#otus_RA_filt_SHD = otus_RA_filt[otus_RA_filt.index.str.contains('D0')]

# LOG TRANSFORMATION OF OTU TAB RELATIVE ABUNDANCES

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

log_otus_all_tab, otus_all_rows, otus_all_cols, otus_all_F = otus_transform(otus_RA_filt)

# COMPUTE AND PLOT BETA DIVERSITY WITH SKBIO
# Simplify metadata for plotting
variable = 'env_classification'
mask_unique_var = metadata_filt[variable].value_counts() <= 2
unique_var = mask_unique_var[mask_unique_var].index
metadata_filt[variable] = metadata_filt[variable].apply(lambda x: 'n<=2' if x in unique_var else x)

# Compute beta diversity
bc_div = skbio.diversity.beta_diversity('braycurtis', log_otus_all_tab, ids=otus_all_rows, validate=True)
bc_pcoa = skbio.stats.ordination.pcoa(bc_div)

# 3D plot for quick visualization
bc_pcoa.plot(metadata_filt,variable,axis_labels=('PC 1', 'PC 2', 'PC 3'),
             title=variable, cmap='tab20', s=24)
plt.tight_layout()
plt.show()

#PLOT PCOA ENV_CLASSIFICATION WITH MATPLOTLIB

df = bc_pcoa.samples  # extract the df
pcoa = pd.merge(df, metadata_filt['env_classification'], left_index=True, right_index=True)
pcoa['env_classification'] = pcoa['env_classification'].replace('Dog others', 'Dogs undet')

# Create a color palette
order = ['Dog Pet', 'Dog Colony', 'Dog Shelter', 'Dog Free_roaming', 'Dogs undet', 'Dog Ancient', 'Wild Canid']
palette = sns.color_palette("Dark2", len(order))
palette[order.index('Dogs undet')] = 'grey'
palette[order.index('Wild Canid')] = '#e6ab02'
palette[order.index('Dog Ancient')] = '#a6761d'
col_dict = dict(zip(order, palette))

# Plot 2D scatterplot
width_mm = 130
height_mm = 60
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)

for idx, row in pcoa.iterrows():
    if 'D0' in idx and row['env_classification'] == 'Dog Pet':
        marker = '^'  # Use triangle marker for SHD cohort
    else:
        marker = 'o'  # Use circle marker for other datasets
        ax.scatter(row['PC1'], row['PC2'], c=col_dict[row['env_classification']], s=20, alpha=0.4, marker=marker)

# SHD plotted last
SHD_data = pcoa[(pcoa['env_classification'] == 'Dog Pet') & (pcoa.index.str.contains('D0'))]
ax.scatter(SHD_data['PC1'], SHD_data['PC2'], c=col_dict['Dog Pet'], s=20, alpha=0.4, marker='^')

# Name axis
ax.set_xlabel('PCo 1', fontsize=10)
ax.set_ylabel('PCo 2', fontsize=10)

# Create legend
handles = [Patch(facecolor=col_dict[name]) for name in order]
plt.legend(handles, order, title=None,
           bbox_to_anchor=(0.61, 0.78), bbox_transform=plt.gcf().transFigure)

plt.tick_params(labelleft=False, labelbottom=False, labeltop=False)
plt.subplots_adjust(left=0.1,right=0.6)
ax.tick_params(which='both', bottom=False, left=False, top=False, right=False)
sns.despine()

# Show/Save figure
#plt.show()
pcoa_path = ('intermediate-outputs/figures/PCoA_all_env_classification_2D_'+ otu_marker +'.svg')
plt.savefig(pcoa_path)

# PLOT PCOA OTHER VARIABLES WITH MATPLOTLIB
variable = 'Animal_age_simplified'
df = bc_pcoa.samples # extract the df
pcoa = pd.merge(df,metadata_filt[variable],left_index=True,right_index=True)

# Create a color palette
unique_var = pcoa[variable].unique()
palette = sns.color_palette("tab20", len(unique_var))
col_dict = dict(zip(unique_var, palette))

# Plot 2D scatterplot
fig, ax = plt.subplots(figsize=(7.5,3.5))
scatter = ax.scatter(pcoa['PC1'], pcoa['PC2'],c=pcoa[variable].map(col_dict),s=40,alpha=0.7)
ax.set_xlabel('PCo 1')
ax.set_ylabel('PCo 2')

# Create legend
handles = [Patch(facecolor=col_dict[variable], label=variable) for variable in unique_var]
plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left', title=None)

# Format axis
plt.tick_params(labelleft=False, labelbottom=False, labeltop=False)
ax.tick_params(which='both', bottom=False, left=False, top=False, right=False)
sns.despine()

# Plot figure
plt.tight_layout()
plt.show()