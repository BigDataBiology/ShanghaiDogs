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
from matplotlib.patches import Patch
import seaborn as sns

os.chdir('/data/Projects/ShanghaiDogs')

# IMPORT INPUT FILES
metadata = pd.read_csv('data/ShanghaiDogsMetadata/ALL_canid_felid_metadata.csv', \
                       sep=',', encoding= 'unicode_escape', index_col=0)
metadata_rep = metadata.query('Genus == "Canis" and Representative=="Yes"')

#otu_tab = pd.read_csv('intermediate-outputs/singlem_profiling/beta-div/all-otu-table.S3.5.rib_prot_S2_rpsB.ebd', \
#                      sep='\t', index_col=0)
#otu_tab = pd.read_csv('intermediate-outputs/singlem_profiling/beta-div/all-otu-table.S3.1.ribosomal_protein_L2_rplB.ebd', \
#                      sep='\t', index_col=0)
#otu_tab = pd.read_csv('intermediate-outputs/singlem_profiling/beta-div/all-otu-table.S3.54.serS.ebd', \
#                      sep='\t', index_col=0)
otu_tab = pd.read_csv('intermediate-outputs/singlem_profiling/beta-div/all-otu-table.S3.18.EIF_2_alpha.ebd', \
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
otus_all_RA.to_csv('intermediate-outputs/singlem_profiling/otus_tab/ALL_OTU_RA_filt_S3.5.rib_prot_S2_rpsB.csv')

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

# Filter metadata to representative samples & remove low mean otus
samples_id = list(otus_all_RA_filt.columns) # Removed samples with low mean otus
metadata_filt = metadata_rep[metadata_rep.index.isin(samples_id)]
samples_id = list(metadata_filt.index)

# Filter otu_tab for Maaslin (counts! before RA calculation)
otus_count_filt = otus_all_filt_wo0[otus_all_filt_wo0.index.isin(samples_id)]
otus_to_keep = list(otus_all_RA_filt.index) # after removal of low mean abundance OTUs
otus_count_filt = otus_count_filt.T
otus_count_filt_abd = otus_count_filt[otus_count_filt.index.isin(otus_to_keep)]
# Save files for Maaslin2
otus_count_filt_abd.to_csv('intermediate-outputs/singlem_profiling/otus_tab/REP_OTU_count_filt_S3.5.rib_prot_S2_rpsB.csv')
metadata_filt.to_csv('intermediate-outputs/singlem_profiling/otus_tab/metadata_filt.csv')

# Filter metadata and RA otu_tab to contain info for only the samples that will be in the PCoA plot
otus_RA_filt = otus_all_RA_filt.T
otus_RA_filt_2 = otus_RA_filt[otus_RA_filt.index.isin(samples_id)]
otus_RA_filt_SHD = otus_RA_filt[otus_RA_filt.index.str.contains('D0')]

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

log_otus_all_tab, otus_all_rows, otus_all_cols, otus_all_F = otus_transform(otus_RA_filt_2)
log_otus_SHD_tab, otus_SHD_rows, otus_SHD_cols, otus_SHD_F = otus_transform(otus_RA_filt_SHD)

### COMPUTE AND PLOT BETA DIVERSITY WITH SKBIO

# Simplify metadata for plotting
variable = 'Breed_subspecies'
mask_unique_var = metadata_filt[variable].value_counts() <= 2
unique_var = mask_unique_var[mask_unique_var].index
metadata_filt[variable] = metadata_filt[variable].apply(lambda x: 'n<=2' if x in unique_var else x)

# Compute beta diversity
bc_div = skbio.diversity.beta_diversity('braycurtis', log_otus_all_tab, ids=otus_all_rows, validate=True)
bc_pcoa = skbio.stats.ordination.pcoa(bc_div)
bc_pcoa.plot(metadata_filt,'Breed_subspecies',axis_labels=('PC 1', 'PC 2', 'PC 3'),
             title='Samples colored by Breed', cmap='tab20', s=24)
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

### PCOA PLOTS WITH MATPLOTLIB
a = bc_pcoa.samples  # extract the df
a2 = pd.merge(a, metadata_rep['env_classification'], left_index=True, right_index=True)

# Updated order with 'Dogs unknown' instead of 'Dog others'
order = ['Dog Pet', 'Dog Colony', 'Dog Shelter', 'Dog Street', 'Wild Canid', 'Dog Ancient','Dogs unknown']
palette = sns.color_palette("Dark2", len(order))

# Replace the color for 'Dogs unknown' with grey
palette[order.index('Dogs unknown')] = 'grey'
palette[order.index('Dog Street')] = '#EF8DB0'
palette[order.index('Wild Canid')] = '#1f78b4'
palette[order.index('Dog Ancient')] = '#e6ab02'

col_dict = dict(zip(order, palette))

# Update the 'env_classification' column in a2 to match the col_dict keys
a2['env_classification'] = a2['env_classification'].replace('Dog others', 'Dogs unknown')

# Plot 2D scatterplot with larger and slightly transparent markers
fig, ax = plt.subplots(figsize=(6, 3.5))

# Plot each point individually with appropriate marker and size
for idx, row in a2.iterrows():
    if 'D0' in idx and row['env_classification'] == 'Dog Pet':
        marker = '^'  # Use triangle marker for entries with 'D0' in index and 'Dog Pet' classification
    else:
        marker = 'o'  # Use circle marker for other entries
    ax.scatter(row['PC1'], row['PC2'], c=col_dict[row['env_classification']], s=40, alpha=0.6, marker=marker)

# Ensure 'Pet dogs' are plotted last (in front of other points)
pet_dogs_data = a2[(a2['env_classification'] == 'Dog Pet') & (a2.index.str.contains('D0'))]
ax.scatter(pet_dogs_data['PC1'], pet_dogs_data['PC2'], c=col_dict['Dog Pet'], s=40, alpha=0.3, marker='^')

ax.set_xlabel('PCo 1')
ax.set_ylabel('PCo 2')

# Create legend with "Dogs unknown" last
handles = [Patch(facecolor=col_dict[name]) for name in order]
plt.legend(handles, order, title=None,
           bbox_to_anchor=(0.95, 0.7), bbox_transform=plt.gcf().transFigure)
plt.tick_params(labelleft=False, labelbottom=False, labeltop=False)
plt.subplots_adjust(right=0.7)

ax.tick_params(which='both', bottom=False, left=False, top=False, right=False)
sns.despine()
#plt.show()
#plt.savefig('analysis/figures/PCoA_all_env_classification_2D_S3.5.rib_prot_S2_rpsB.svg')
#plt.savefig('analysis/figures/PCoA_all_env_classification_2D_S3.1.ribosomal_protein_L2_rplB.svg')
#plt.savefig('analysis/figures/PCoA_all_env_classification_2D_S3.54.serS.svg')
plt.savefig('analysis/figures/PCoA_all_env_classification_2D_S3.14.hisS.svg')

### PCOA PLOTS WITH MATPLOTLIB STUDY

a = bc_pcoa.samples # extract the df
a2 = pd.merge(a,metadata_rep['Study'],left_index=True,right_index=True)

# Create a color palette
unique_studies = a2['Study'].unique()
palette = sns.color_palette("tab20", len(unique_studies))
col_dict = dict(zip(unique_studies, palette))

#Plot 2D scatterplot
fig, ax = plt.subplots(figsize=(7.5,3.5))
scatter = ax.scatter(a2['PC1'], a2['PC2'],c=a2['Study'].map(col_dict),s=40,alpha=0.7)
ax.set_xlabel('PCo 1')
ax.set_ylabel('PCo 2')

# Create legend handles
handles = [Patch(facecolor=col_dict[study], label=study) for study in unique_studies]
plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left', title=None)

plt.tick_params(labelleft=False, labelbottom=False, labeltop=False)
plt.subplots_adjust(right=0.8)

ax.tick_params(which='both', bottom=False, left=False, top=False, right=False)
sns.despine()
plt.tight_layout()
plt.show()