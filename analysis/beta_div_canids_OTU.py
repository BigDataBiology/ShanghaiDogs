# -*- coding: utf-8 -*-
"""
Created on Mon 08 Apr 2024 16:16:51
@author: Anna Cusco
"""

import pandas as pd
import numpy as np
import skbio
from skbio.stats.distance import DistanceMatrix, permanova, permdisp
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import csv

plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

# IMPORT INPUT FILES

# Metadata
metadata = pd.read_csv('data/ShanghaiDogsMetadata/REP_canid_metadata.csv', \
                       sep=',', encoding= 'unicode_escape', index_col=0)
metadata['env_classification'] = metadata['env_classification'].str.replace(' captive','')

# 2) SRA run to biosample file
SRA_metadata = pd.read_csv('external-data/data/dog_microbiome_archive_otu_tables/SRR_Acc_List_Metadata.txt', \
                               sep='|',skiprows=[1],index_col=0)
run_to_biosample = SRA_metadata[[' biosample      ']].reset_index()
run_to_biosample.columns = ['run','biosample']
run_to_biosample['run'] = run_to_biosample['run'].str.replace(' ', '', regex=True)
run_to_biosample['biosample'] = run_to_biosample['biosample'].str.replace(' ', '', regex=True)
run_to_biosample = run_to_biosample.set_index('run')

# 3) SingleM marker genes
singlem_markers = pd.read_csv('intermediate-outputs/singlem_profiling/beta-div/singlem_metapackage_describe.txt', \
                              sep='\t', encoding= 'unicode_escape', index_col=0)
singlem_markers_arch_bact = singlem_markers[singlem_markers['target_domains']=='Bacteria,Archaea']
otu_marker_ls = singlem_markers_arch_bact.index.to_list()

# FORMATTING & FILTERING OTU TABLE - Function definition
# Adapt index names: rename run_ids to Biosample_ids
def reformat_idx_otu_tab(otu_tab):
    """
    reformat indexes from otu_tab in order to later match metadata
    """
    # Reformat index names: remove suffixes
    otu_tab.index = otu_tab.index.str.replace('_350', '')
    otu_tab.index = otu_tab.index.str.replace('_1', '')
    otu_tab.index = otu_tab.index.str.replace('_EKDN.*', '', regex=True)
    for idx in otu_tab.index:
        if (idx.startswith('SRR') or idx.startswith('ERR')) and idx in run_to_biosample.index:
            otu_tab.rename(index={idx: run_to_biosample.loc[idx, 'biosample']}, inplace=True)
    return(otu_tab)

# 1) Define filtering functions
# Filter out OTUs with low total counts and low prevalence
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

# Filter out SAMPLES with low OTU count
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

# Remove anything that is 0 on rows or columns for final filtering
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

# Calculate relative abundance of OTU table
def rel_ab_otu(otu_tab):
    """
    check OTU table orientation, tranpose when necessary,
    and calculate relative abundance of OTUs/sample.
    Index (rows) should contain the samples
    """
    otu_tab['Sum'] = otu_tab.sum(axis=1)
    otus_RA = otu_tab.iloc[:, :].div(otu_tab.Sum, axis=0)
    otus_RA.drop(['Sum'], axis=1, inplace=True)
    otu_tab.drop(['Sum'], axis=1, inplace=True)
    return otus_RA

# Remove OTUs with low mean relative abundance
def rm_low_mean_otus(OTUs_tab,min):
    """
    transpose OTU table (OTU_ids should be on the rows)
    calculate mean mOTU relative abundance and remove those
    motus with a lower mean abundance (min).
    """
    OTUs_tab = OTUs_tab.T
    print(OTUs_tab.index[0] + ': this should be an OTU ID')
    OTUs_tab['Mean'] = OTUs_tab.mean(axis=1)
    OTUs_tab_filt = OTUs_tab.loc[OTUs_tab['Mean'].between(min, 1)]
    OTUs_tab_filt = OTUs_tab_filt.drop(['Mean'], axis=1)
    return OTUs_tab_filt

# log transformation of otu_tab relative abundances
def otus_transform(OTUs_tab):
    """
    Multiply relative abundances by 10^6. Calculate a
    pseudocount by taking the minimum value from OTU table
    and divide it by 2. Substitute zeros for the pseudocount.
    Perform log10 transformation of the matrix.

    Required in OTU_tab: OTUs IDs need to be on the columns.
    """
    OTUs_tab = OTUs_tab.T
    print(OTUs_tab.index[0] + ': this should be an OTU ID')
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

# DEFINING COMMON SAMPLE_LIST (across all OTU tables)
common_samples = set(metadata.index)

for m in otu_marker_ls:
    print('Updating number of samples... ')
    print(len(common_samples))
    # 1) import otu_table
    otu_marker = m
    input_path = ('intermediate-outputs/singlem_profiling/beta-div/unifrac-otu/all-otu-table.'+otu_marker+'.ebd')
    otu_tab = pd.read_csv(input_path, sep='\t', index_col=0)

    # 2) Filter the otu table
    otu_rf = reformat_idx_otu_tab(otu_tab)
    otus_filt_0, otus_filt_rm_otus = filt_low_OTU(otu_rf, 10, 0.02) # filt for total OTU count and prev
    otus_filt_1, otus_filt_rm_samples = filt_samples_low_OTU(otus_filt_0, 200) # filt for min sum OTU / sample
    otus_filt, rm_all_items = remove_0_sum(otus_filt_1)  # remove 0s
    if 'D024' in otus_filt.index:
        otus_filt.drop('D024', axis=0, inplace=True)  # in ATBs treatment

    # 3) Calculate RA and log transform
    otus_RA = rel_ab_otu(otus_filt)
    otus_RA_filt = rm_low_mean_otus(otus_RA, 0.00005) # optional, we could include all species min=0
    samples_id = set(otus_RA_filt.columns)

    # 4) Adjust samples from metadata and otu_tab
    common_samples &= samples_id

# Now we can proceed to compute beta div and calculate the median distances
# (considering the same exact comparisons across samples)

# FORMATTING, FILTERING, & COMPUTE BETA DIV
print('Final number of samples for all included marker genes is:')
print(len(common_samples))

ls_samples = list(common_samples)
with open('intermediate-outputs/singlem_profiling/beta-div/samples_ls.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(ls_samples)

metadata_filt = metadata[metadata.index.isin(ls_samples)]
distance_matrices = []

for m in otu_marker_ls:
    # 1) import otu_table
    otu_marker = m
    input_path = ('intermediate-outputs/singlem_profiling/beta-div/unifrac-otu/all-otu-table.'+otu_marker+'.ebd')
    otu_tab = pd.read_csv(input_path, sep='\t', index_col=0)

    # 2) Filter the otu table
    otu_rf = reformat_idx_otu_tab(otu_tab)
    otus_filt_0, otus_filt_rm_otus = filt_low_OTU(otu_rf, 10, 0.02) # filt for total OTU count and prev
    otus_filt_1, otus_filt_rm_samples = filt_samples_low_OTU(otus_filt_0, 200) # filt for min sum OTU / sample
    otus_filt, rm_all_items = remove_0_sum(otus_filt_1)  # remove 0s
    if 'D024' in otus_filt.index:
        otus_filt.drop('D024', axis=0, inplace=True)  # in ATBs treatment

    # 3) Calculate RA and log transform
    otus_RA = rel_ab_otu(otus_filt)
    otus_RA_filt = rm_low_mean_otus(otus_RA, 0.00005) # optional, we could include all species min=0 or filter minimally 0.00005

    # 4) Adjust samples from metadata and otu_tab
    otus_RA_out = otus_RA_filt.T
    otus_RA_out = otus_RA_out[otus_RA_out.index.isin(ls_samples)]
    print(len(otus_RA_out))

    # 5) Log transformation of filtered otus_RA table
    log_otus_tab, samples_ids_rows, otus_ids_cols, otus_tab_F = otus_transform(otus_RA_out)

    # 6) Compute beta_diversity
    bc_div = skbio.diversity.beta_diversity('braycurtis', log_otus_tab, ids=samples_ids_rows, validate=True)

    # 7) Store all distance matrices, and compute the median
    distance_matrices.append(bc_div)

# Compute the median distance matrix across all marker genes
distance_matrices_array = np.array([dm.data for dm in distance_matrices])  # Extract distance matrix data
median_distance_array = np.median(distance_matrices_array, axis=0)
median_distance_matrix = DistanceMatrix(median_distance_array, ids=distance_matrices[0].ids)

# COMPUTE AND PLOT BETA DIVERSITY WITH SKBIO
# Simplify metadata for plotting
variable = 'env_classification'
mask_unique_var = metadata_filt[variable].value_counts() <= 2
unique_var = mask_unique_var[mask_unique_var].index
metadata_filt[variable] = metadata_filt[variable].apply(lambda x: 'n<=2' if x in unique_var else x)

# Compute beta diversity
bc_pcoa = skbio.stats.ordination.pcoa(median_distance_matrix)

# 3D plot for quick visualization
bc_pcoa.plot(metadata_filt,variable,axis_labels=('PC 1', 'PC 2', 'PC 3'),
             title=variable, cmap='tab20', s=24)
plt.tight_layout()
plt.show()

#PLOT PCOA ENV_CLASSIFICATION WITH MATPLOTLIB

df = bc_pcoa.samples  # extract the df
pcoa = pd.merge(df, metadata_filt['env_classification'], left_index=True, right_index=True)
pcoa['env_classification'] = pcoa['env_classification'].replace('Dog others', 'Dog undet')

# Create a color palette
order = ['Dog Pet', 'Dog Colony', 'Dog Shelter', 'Dog Free_roaming', 'Wild Canid', 'Dog undet']
palette = sns.color_palette("Dark2", len(order))
palette[order.index('Dog Shelter')] = '#e6ab02'
palette[order.index('Dog Free_roaming')] = '#7570b3'
palette[order.index('Dog undet')] = 'grey'
palette[order.index('Wild Canid')] = '#e7298a'

#palette[order.index('Dog Ancient')] = '#a6761d'
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
plt.savefig('intermediate-outputs/figures/PCoA-2D-CANID-env_class-median_markers-filt.svg')


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

### Check grouping statistics

print(median_distance_matrix)
ls_variables = ['env_classification', 'Animal_age_simplified',
                      'Sex', 'Study', 'Size_class']

# Define filtering function
def filter_distance_matrix(distance_matrix, metadata, variable):
    """
    Filters a distance matrix by excluding samples based on
    conditions applied to metadata.
    """
    # Study list to keep
    study_ls = ['This_study', 'Coelho_2018', 'Allaway_2020_dogs', 'Bai_2023', 'Wang_2019_dogs', 'Berlin_cohort']

    # Filter metadata
    filtered_metadata = metadata.loc[~metadata[variable].isin
    (['Dog others', 'Dog Shelter', 'unclass', 'unknown', 'Unknown'])]
    if variable in ['Animal_age_simplified', 'Sex', 'Size_class']:
        filtered_metadata = filtered_metadata[filtered_metadata['env_classification'] != 'Wild Canid']
    if variable == 'Study':
        filtered_metadata = filtered_metadata[filtered_metadata[variable].isin(study_ls)]
    print(filtered_metadata[variable].unique())

    # Get filtered sample IDs and filter distance matrix
    filtered_ids = filtered_metadata.index.intersection(distance_matrix.ids)
    filtered_distance_matrix = distance_matrix.filter(filtered_ids, strict=False)
    return filtered_distance_matrix, filtered_metadata

# Initialize results_dataframe

permdisp_res = []
permanova_res = []

for variable in ls_variables:
    filt_distance_matrix, filt_metadata = filter_distance_matrix(median_distance_matrix, metadata_filt, variable)
    column = variable
    metadata_df = filt_metadata[[variable]]
    print("Running PERMDISP for " + variable)
    permdisp_result = permdisp(filt_distance_matrix, metadata_df[variable], permutations=999)

    # Store PERMDISP results in a dictionary
    permdisp_res.append({
        "Variable": variable,
        "Method": permdisp_result["method name"],
        "Test Statistic Name": permdisp_result["test statistic name"],
        "Test Statistic": permdisp_result["test statistic"],
        "p-value": permdisp_result["p-value"],
        "Sample Size": permdisp_result["sample size"],
        "Number of Groups": permdisp_result["number of groups"],
        "Number of Permutations": permdisp_result["number of permutations"]})

    print("Running PERMANOVA for " + variable)
    permanova_result = permanova(filt_distance_matrix, metadata_df[variable], permutations=999)

    # Store PERMANOVA results in a dictionary
    permanova_res.append({
        "Variable": variable,
        "Method": permanova_result["method name"],
        "Test Statistic Name": permanova_result["test statistic name"],
        "Test Statistic": permanova_result["test statistic"],
        "p-value": permanova_result["p-value"],
        "Sample Size": permanova_result["sample size"],
        "Number of Groups": permanova_result["number of groups"],
        "Number of Permutations": permanova_result["number of permutations"]})

# Convert results to DataFrames for easier visualization or saving
permdisp_df = pd.DataFrame(permdisp_res)
permanova_df = pd.DataFrame(permanova_res)

# Multiple test correction
permdisp_df["Adjusted p-value"] = multipletests(permdisp_df["p-value"], method="fdr_bh")[1]
permanova_df["Adjusted p-value"] = multipletests(permanova_df["p-value"], method="fdr_bh")[1]

permdisp_df.to_csv('intermediate-outputs/singlem_profiling/beta-div/permdisp_res_global.csv')
permanova_df.to_csv('intermediate-outputs/singlem_profiling/beta-div/permanova_res_global.csv')

## Compute pairwise comparisons using PERMANOVA
pairwise_all = []

for variable in ls_variables:
    filt_distance_matrix, filt_metadata = filter_distance_matrix(median_distance_matrix, metadata_filt, variable)
    column = variable
    metadata_df = filt_metadata[[variable]]
    # Pairwise Comparisons (optional, for multiple groups)
    unique_groups = metadata_df[variable].unique()
    pairwise_results = []

    for i, group1 in enumerate(unique_groups):
        for group2 in unique_groups[i + 1:]:
            print(f"Running pairwise PERMANOVA for {group1} vs {group2}...")
            pairwise_mask = metadata_df[variable].isin([group1, group2])
            pairwise_ids = metadata_df.index[pairwise_mask]
            pairwise_dm = filt_distance_matrix.filter(pairwise_ids)
            pairwise_gv = metadata_df[variable][pairwise_mask]
            pairwise_result = permanova(pairwise_dm, pairwise_gv, permutations=999)

            pairwise_results.append({
                "Variable": variable,
                "Group 1": group1,
                "Num samples group 1": len(metadata_df[metadata_df[variable]==group1]),
                "Group 2": group2,
                "Num samples group 2": len(metadata_df[metadata_df[variable] == group2]),
                "Test Statistic": pairwise_result["test statistic"],
                "p-value": pairwise_result["p-value"],
                "Sample Size": pairwise_result["sample size"],
                "Number of Groups": pairwise_result["number of groups"],
                "Number of Permutations": pairwise_result["number of permutations"]})

    # Convert results to a DataFrame and perform multiple test corr (for the current variable)
    pairwise_results_df = pd.DataFrame(pairwise_results)
    pairwise_results_df["Adjusted p-value"] = multipletests(pairwise_results_df["p-value"], method="fdr_bh")[1]
    pairwise_all.append(pairwise_results_df)

# Merge all pairwise comparison DataFrames into a single DataFrame
combined_pairwise_df = pd.concat(pairwise_all, ignore_index=True)
combined_pairwise_df['Adjusted p-value (global)'] = multipletests(combined_pairwise_df['p-value'], method="fdr_bh")[1]
combined_pairwise_df.to_csv('intermediate-outputs/singlem_profiling/beta-div/permanova_pairwise_results.csv'
                            , sep=',', index=False)


## Compute pairwise comparisons using PERMDISP
# (to assess between which groups there is variability in dispersion rather than centroids)

permdisp_pw_all = []

for variable in ls_variables:
    filt_distance_matrix, filt_metadata = filter_distance_matrix(median_distance_matrix, metadata_filt, variable)
    column = variable
    metadata_df = filt_metadata[[variable]]
    # Pairwise Comparisons (optional, for multiple groups)
    unique_groups = metadata_df[variable].unique()
    pairwise_results = []

    for i, group1 in enumerate(unique_groups):
        for group2 in unique_groups[i + 1:]:
            print(f"Running pairwise PERMDISP for {group1} vs {group2}...")
            pairwise_mask = metadata_df[variable].isin([group1, group2])
            pairwise_ids = metadata_df.index[pairwise_mask]
            pairwise_dm = filt_distance_matrix.filter(pairwise_ids)
            pairwise_gv = metadata_df[variable][pairwise_mask]
            pairwise_result = permdisp(pairwise_dm, pairwise_gv, permutations=999)

            pairwise_results.append({
                "Variable": variable,
                "Group 1": group1,
                "Num samples group 1": len(metadata_df[metadata_df[variable]==group1]),
                "Group 2": group2,
                "Num samples group 2": len(metadata_df[metadata_df[variable] == group2]),
                "Test Statistic": pairwise_result["test statistic"],
                "p-value": pairwise_result["p-value"],
                "Sample Size": pairwise_result["sample size"],
                "Number of Groups": pairwise_result["number of groups"],
                "Number of Permutations": pairwise_result["number of permutations"]})

    # Convert results to a DataFrame and perform multiple test corr (for the current variable)
    pairwise_results_df = pd.DataFrame(pairwise_results)
    pairwise_results_df["Adjusted p-value"] = multipletests(pairwise_results_df["p-value"], method="fdr_bh")[1]
    permdisp_pw_all.append(pairwise_results_df)

# Merge all pairwise comparison DataFrames into a single DataFrame
combined_permdisp_pw_df = pd.concat(permdisp_pw_all, ignore_index=True)
combined_permdisp_pw_df['Adjusted p-value (global)'] = multipletests(combined_pairwise_df['p-value'], method="fdr_bh")[1]
combined_permdisp_pw_df.to_csv('intermediate-outputs/singlem_profiling/beta-div/permdisp_pairwise_results.csv'
                            , sep=',', index=False)

### PERMANOVA & PERMDISP results

global_permdisp = pd.read_csv('intermediate-outputs/singlem_profiling/beta-div/permdisp_res_global.csv',
                              sep=',', index_col=0)
global_permanova = pd.read_csv('intermediate-outputs/singlem_profiling/beta-div/permanova_res_global.csv',
                               sep=',', index_col=0)

pw_permdisp = pd.read_csv('intermediate-outputs/singlem_profiling/beta-div/permdisp_pairwise_results.csv',
                                      sep=',')
pw_permanova = pd.read_csv('intermediate-outputs/singlem_profiling/beta-div/permanova_pairwise_results.csv',
                                      sep=',')

# PLOT PCOAs FOR FIGURE (env_classification)

filt_distance_matrix, filt_metadata = filter_distance_matrix(median_distance_matrix,
                                                             metadata_filt,
                                                             'env_classification')
bc_pcoa = skbio.stats.ordination.pcoa(filt_distance_matrix)
eigenvalues = bc_pcoa.eigvals
explained_variance = (eigenvalues / eigenvalues.sum()) * 100

# Prepare the dataframe to plot
df = bc_pcoa.samples  # extract the df
pcoa = pd.merge(df, metadata_filt['env_classification'], left_index=True, right_index=True)
pcoa['env_classification'] = pcoa['env_classification'].replace('Dog others', 'Dog undet')

# Plot 2D scatterplot
width_mm = 130
height_mm = 60
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)

handles = []
labels = []
order = ['Dog Pet (This study)', 'Dog Pet', 'Dog Colony', 'Dog Free_roaming', 'Wild Canid']

for idx, row in pcoa.iterrows():
    if 'D0' in idx and row['env_classification'] == 'Dog Pet':
        marker = '^'  # Use triangle marker for SHD cohort
        label = 'Dog Pet (This study)'
    else:
        marker = 'o'  # Use circle marker for other datasets
        label = row['env_classification']
        ax.scatter(row['PC1'], row['PC2'], c=col_dict[row['env_classification']], s=20, alpha=0.4, marker=marker)

    # Add the scatter plot for the legend
    if label not in labels:  # Avoid adding duplicate labels
        handles.append(ax.scatter([], [], c=col_dict.get(label, 'grey'), marker=marker, s=20, alpha=0.4))
        labels.append(label)

# SHD plotted last
SHD_data = pcoa[(pcoa['env_classification'] == 'Dog Pet') & (pcoa.index.str.contains('D0'))]
ax.scatter(SHD_data['PC1'], SHD_data['PC2'], c=col_dict['Dog Pet'], s=22, alpha=0.4, marker='^')

# Name axis
xlab = f'Axis 1 ({explained_variance[0]:.1f}%)'
ylab = f'Axis 2 ({explained_variance[1]:.1f}%)'
ax.set_xlabel(xlab, fontsize=10)
ax.set_ylabel(ylab, fontsize=10)
ax.tick_params(which='both', bottom=False, left=False, top=False, right=False)
plt.tick_params(labelleft=False, labelbottom=False, labeltop=False)

# Create legend
handles = [h for l, h in sorted(zip(labels, handles), key=lambda x: order.index(x[0])) if l in order]
labels = [l for l in order]
plt.legend(handles=handles, labels=labels, title="",
           bbox_to_anchor=(0.61, 0.78), bbox_transform=plt.gcf().transFigure)
plt.subplots_adjust(left=0.1,right=0.6)


# Show/Save figure
sns.despine()
plt.tight_layout()
#plt.show()
plt.savefig('intermediate-outputs/figures/PCoA-2D-CANID-env_class_4cat-median_markers-filt.svg')
