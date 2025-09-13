# -*- coding: utf-8 -*-
"""
Created on Thu 21 Nov 2024 12:21:06
@author: Anna Cusco
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from skbio.diversity import alpha_diversity
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from itertools import combinations

plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

# IMPORT INPUT FILES

# 1) Metadata
metadata = pd.read_csv('data/ShanghaiDogsMetadata/REP_canid_metadata.csv', \
                       sep=',', encoding= 'unicode_escape', index_col=0)
metadata['env_classification'] = metadata['env_classification'].str.replace(' captive','')
metadata['Sex'] = metadata['Sex'].str.replace('Female','female')
metadata['Sex'] = metadata['Sex'].str.replace('Male','male')
metadata['Study'] = metadata['Study'].str.replace(r'_dogs|_Canidae|_coprolite|_global_dog', '',regex=True)

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

# 4) Initialize results df
summary_shannon = pd.DataFrame()
summary_simpson = pd.DataFrame()
kruskal_all = []
mannwhitneyu_all = []

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

# FORMATTING, FILTERING, & COMPUTE ALPHA DIV

for m in otu_marker_ls:
    # 1) import otu_table
    otu_marker = m
    input_path = ('intermediate-outputs/singlem_profiling/beta-div/unifrac-otu/all-otu-table.'+otu_marker+'.ebd')
    otu_tab = pd.read_csv(input_path, sep='\t', index_col=0)

    # 2) Filter the otu table
    otu_rf = reformat_idx_otu_tab(otu_tab)
    otus_filt_0, otus_filt_rm_otus = filt_low_OTU(otu_rf, 10, 0.02) # filt for total OTU count and prev
    otus_filt_1, otus_filt_rm_samples = filt_samples_low_OTU(otus_filt_0, 200) # filt for min sum OTU / sample
    otus_filt, rm_all_items = remove_0_sum(otus_filt_1) # remove 0s
    if 'D024' in otus_filt.index:
        otus_filt.drop('D024', axis=0, inplace=True)  # in ATBs treatment

    # Calculate RA
    otus_filt_RA = rel_ab_otu(otus_filt)

    # Calculate shannon and simpson indexes
    shannon = alpha_diversity('shannon', otus_filt_RA, otus_filt_RA.index).to_frame(name='shannon_' + otu_marker)
    simpson = alpha_diversity('simpson', otus_filt_RA, otus_filt_RA.index).to_frame(name='simpson_' + otu_marker)

    summary_shannon = pd.merge(summary_shannon, shannon, left_index=True, right_index=True, how='outer')
    summary_simpson = pd.merge(summary_simpson, simpson, left_index=True, right_index=True, how='outer')

# Keep only representative samples
summary_shannon = summary_shannon[summary_shannon.index.isin(metadata.index)]
summary_simpson = summary_simpson[summary_simpson.index.isin(metadata.index)]

# Remove samples that alpha div is nan for more than half of the chosen marker genes (>6)
nan_counts_shannon = summary_shannon.isna().sum(axis=1)
nan_counts_simpson = summary_simpson.isna().sum(axis=1)

filtered_shannon = summary_shannon[nan_counts_shannon <= 6]
filtered_simpson = summary_simpson[nan_counts_simpson <= 6]

# Calculate the median alpha diversity values for all archaeal and bacterial marker genes
filtered_shannon['median'] = filtered_shannon.median(axis=1)
filtered_simpson['median'] = filtered_simpson.median(axis=1)

# Merge with metadata variables that want to assess significance
shannon_w_metadata = pd.merge(filtered_shannon['median'],metadata[['Size_class', 'Animal_age_simplified',
                                                                  'Sex', 'Study', 'env_classification']],
                              left_index=True, right_index=True)

simpson_w_metadata = pd.merge(filtered_simpson['median'],metadata[['Size_class', 'Animal_age_simplified',
                                                                  'Sex', 'Study', 'env_classification']],
                              left_index=True, right_index=True)

# Calculate statistical significance (KW + MANN WHITNEY U) + multiple-testing corrections

# Evaluated variables
variable_ls = ['Size_class', 'Animal_age_simplified',
               'Sex', 'Study', 'env_classification']

# Filter data function
def filter_data(df, variable):
    study_ls = ['This_study', 'Coelho_2018', 'Allaway_2020', 'Bai_2023', 'Wang_2019']
    filtered = df.loc[~df[variable].isin(['Dog others', 'Dog Shelter', 'unclass', 'unknown', 'Unknown'])]
    if variable in ['Animal_age_simplified', 'Sex', 'Size_class']:
        filtered = filtered[filtered['env_classification'] != 'Wild Canid']
    if variable == 'Study':
        filtered = filtered[filtered[variable].isin(study_ls)]
    return filtered

# Kruskal-Wallis test

shannon_KW = []
simpson_KW = []

for variable in variable_ls:
    # Shannon index
    shannon_filt = filter_data(shannon_w_metadata, variable)
    shannon_grouped = [shannon_filt['median'][shannon_filt[variable] == g].values
                    for g in shannon_filt[variable].unique()]
    if all(len(group) > 1 for group in shannon_grouped):
        kw_stat, kw_pval = kruskal(*shannon_grouped)
        shannon_KW.append((variable, kw_stat, kw_pval))
    # Simpson index
    simpson_filt = filter_data(simpson_w_metadata, variable)
    simpson_grouped = [simpson_filt['median'][simpson_filt[variable] == g].values
        for g in simpson_filt[variable].unique()]
    if all(len(group) > 1 for group in simpson_grouped):
        kw_stat, kw_pval = kruskal(*simpson_grouped)
        simpson_KW.append((variable, kw_stat, kw_pval))

# Apply multiple testing correction for Kruskal-Wallis
kw_results = shannon_KW #simpson_KW

if kw_results:
    variables, kw_stats, kw_pvalues = zip(*kw_results)
    _, kw_corrected_pvalues, _, _ = multipletests(kw_pvalues, method='fdr_bh')
    kruskal_df = pd.DataFrame({
        'Variable': variables,
        'KW Statistic': kw_stats,
        'Original p-value': kw_pvalues,
        'Corrected p-value': kw_corrected_pvalues
    })
    print("\nKruskal-Wallis Results (Corrected):")
    print(kruskal_df)
    kruskal_df.to_csv('intermediate-outputs/singlem_profiling/alpha-div/KW_shannon_results.csv')
else:
    print("No significant Kruskal-Wallis results.")

# When min total sum/sample is >400, size becomes significant variable (due to the outlier sample, probably)
kruskal_df = pd.read_csv('intermediate-outputs/singlem_profiling/alpha-div/KW_shannon_results_200.csv',index_col=0)

# For significant kruskal-wallis results (corrected p-values), plot boxplot
# Only Shannon, for study and env_classification

kruskal_sign = kruskal_df[kruskal_df['Corrected p-value']<0.05]
variable_ls = list(kruskal_sign['Variable'])

for variable in variable_ls:
    shannon_filt = filter_data(shannon_w_metadata, variable)
    categories = sorted(shannon_filt[variable].unique())
    kw_corr_pvalue = kruskal_sign.loc[kruskal_sign['Variable'] == variable, 'Corrected p-value'].iloc[0]

    # Figure size
    width_mm = 120
    height_mm = 80
    figsize_inch = (width_mm / 25.4, height_mm / 25.4)

    # Plot figure
    fig, ax = plt.subplots(figsize=figsize_inch, gridspec_kw={'top': 0.85, 'bottom': 0.40, 'right': 0.60})
    # Plot boxplot + swarmplot
    sns.swarmplot(data=shannon_filt, x=variable, y='median',
                  palette='Dark2', hue=variable, order=categories,
                  alpha=0.7, size=2.5)
    sns.boxplot(data=shannon_filt, x=variable, y='median',
                color='white', width=0.5, order=categories)
    sns.despine()
    ax.text(0.5, 1.2, f'KW corr p-value = {kw_corr_pvalue:.2e}',
            ha='center', va='center', transform=ax.transAxes, fontsize=10)
    ax.set_xticklabels(categories, rotation=45, ha='right')
    ax.set_xlabel('')
    ax.set_ylabel('Shannon index')
    ax.tick_params(axis='y', rotation=0)
    legend = ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.1), fontsize=10)
    out_path = 'intermediate-outputs/figures/alpha_div/ALL/ALL_'+ variable + '_shannon.svg'
    plt.savefig(out_path)
    #plt.show()

## Continue and compute Mann-Whitney U pairwise tests
pairwise_results = []

# Iterate through significant variables
for _, row in kruskal_sign.iterrows():
    variable = row['Variable']

    # Filter the data for the variable
    shannon_filt = filter_data(shannon_w_metadata, variable)

    # Get unique groups for pairwise comparisons
    groups = sorted(shannon_filt[variable].unique())

    # Perform Mann-Whitney U test for each pair of groups
    for g1, g2 in combinations(groups, 2):
        group1 = shannon_filt['median'][shannon_filt[variable] == g1]
        group2 = shannon_filt['median'][shannon_filt[variable] == g2]

        # Ensure both groups have >1 sample
        if len(group1) > 1 and len(group2) > 1:
            stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
            pairwise_results.append((variable, g1, g2, stat, p))

pairwise_results_df = pd.DataFrame(pairwise_results)
pairwise_results_df.columns = ['Variable','g1','g2','stat','Original p-value']


# Store corrected results
corrected_results = []

# Iterate through significant variables:
for variable in pairwise_results_df['Variable'].unique():
    variable_results = pairwise_results_df.loc[(pairwise_results_df['Variable'] == variable)]
    print(variable)
    print(variable_results)

    if not variable_results.empty:  # Ensure there are results to process
        # Extract p-values for this variable
        variable_pvalues = variable_results['Original p-value'].values

        # Apply multiple testing correction (FDR) within this variable
        _, variable_corrected_pvalues, _, _ = multipletests(variable_pvalues, method='fdr_bh')

        # Combine results with corrected p-values
        corrected_results.extend([
            (*row, corr_p) for row, corr_p in zip(variable_results.to_records(index=False), variable_corrected_pvalues)
        ])

# Convert to DataFrame for easy analysis
mannwhitneyu_all_df = pd.DataFrame(corrected_results, columns=['Variable', 'Group 1', 'Group 2',
                                                               'Statistic', 'Original p-value',
                                                               'Corrected p-value'])

### Simplified alpha div plot for main manuscript
variable = 'env_classification'
shannon_filt = filter_data(shannon_w_metadata, variable)
kw_corr_pvalue = kruskal_sign.loc[kruskal_sign['Variable'] == variable, 'Corrected p-value'].iloc[0]

order = ['Dog Pet', 'Dog Colony', 'Dog Free_roaming', 'Wild Canid']

shannon_filt[variable] = pd.Categorical(
    shannon_filt[variable],
    categories=order,
    ordered=True)

shannon_filt = shannon_filt.sort_values(
    by=[variable],
    ascending=[False])

# Figure size
width_mm = 65
height_mm = 50
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)

sns.swarmplot(data=shannon_filt, x=variable, y='median',
              palette='Dark2', hue=variable, alpha = 0.6,
              size=2.5, order=order)
sns.boxplot(data=shannon_filt, x=variable, y='median',
            color='white', order=order, width=0.6)

ax.text(0.5, 1.1, f'KW corr p-value = {kw_corr_pvalue:.2e}',
        ha='center', va='center', transform=ax.transAxes, fontsize=9)
ax.tick_params(bottom=False,labelbottom=False)
ax.set_ylabel('Shannon index')
ax.set_xlabel('')

ax.get_legend().remove()
sns.despine()
plt.tight_layout()
#plt.show()
plt.savefig('intermediate-outputs/figures/alpha_div/ALL/box-swarmplot_'+ variable + '_shannon.svg')

