# -*- coding: utf-8 -*-

"""
Created on Mon Jan 20 12:26:01 2025
@author: Anna Cusco
"""

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy.stats import kruskal,mannwhitneyu
from statsmodels.stats.multitest import multipletests
from itertools import combinations

os.chdir('/data/Projects/ShanghaiDogs/')
plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

### Import data
## load sample list of final beta-div plot
with open('intermediate-outputs/singlem_profiling/beta-div/samples_ls.csv', 'r') as file:
    data = file.read()
samples_ls = data.split(',')

## import metadata
metadata = pd.read_csv('data/ShanghaiDogsMetadata/REP_canid_metadata.csv', \
                       sep=',', encoding= 'unicode_escape', index_col=0)
metadata['env_classification'] = metadata['env_classification'].str.replace(' captive','')

## import species-level tax profile
tax_profile = pd.read_csv('intermediate-outputs/singlem_profiling/tax-profiles/all-dog-tax-profile-species.tsv', \
                          delimiter='\t', header=0, index_col=0)

# Reformat column and index names
tax_profile.columns = tax_profile.columns.str.replace('_350','')
tax_profile.columns = tax_profile.columns.str.replace('_1','')
tax_profile.columns = tax_profile.columns.str.replace('_EKDN.*','',regex=True)
tax_profile.index = tax_profile.index.str.replace('.*s__','',regex=True)

# Reformat column names: rename run_ids to Biosample_ids
SRA_metadata = pd.read_csv('external-data/data/dog_microbiome_archive_otu_tables/SRR_Acc_List_Metadata.txt', \
                               sep='|',skiprows=[1],index_col=0)
run_to_biosample = SRA_metadata[[' biosample      ']].reset_index()
run_to_biosample.columns = ['run','biosample']
run_to_biosample['run'] = run_to_biosample['run'].str.replace(' ', '', regex=True)
run_to_biosample['biosample'] = run_to_biosample['biosample'].str.replace(' ', '', regex=True)
run_to_biosample = run_to_biosample.set_index('run')

for name in tax_profile.columns:
    if (name.startswith('SRR') or name.startswith('ERR')) and name in run_to_biosample.index:
        tax_profile.rename(columns={name: run_to_biosample.loc[name, 'biosample']}, inplace=True)

### Filter tax-profile to include just the relevant samples
tax_profile_filt = tax_profile[tax_profile.columns[tax_profile.columns.isin(samples_ls)]]
tax_profile_filt = tax_profile_filt.T
metadata_filt = metadata[metadata.index.isin(samples_ls)]
tax_w_metadata = pd.merge(tax_profile_filt,metadata_filt[['env_classification','Animal_age_simplified']],left_index=True,right_index=True)

## Significant species list names according to Maaslin2
env_sign_sp = ["Collinsella sp008014645", "Faecalimonas umbilicata", "Sutterella sp905186105",
               "Bifidobacterium globosum", "Phocaeicola plebeius", "Turicibacter sp002311155"]

# Remove "Blautia_A sp019416265", "Fournierella massiliensis", "Blautia sp900556555" - not sign when KW test
age_sign_sp = ["CAJMNU01 sp905214855", "Butyricicoccus pullicaecorum", "Thomasclavelia spiroformis_A",
               "JAGZHZ01 sp018366495", "UMGS1370 sp900551135", "Faecalibacillus intestinalis",
               "Faecalibacterium sp900540455", "Slackia_A piriformis", "Eisenbergiella sp900539715",
               "Dorea_B phocaeensis"]

### Env_classification: Plot boxplots and add pairwise comparisons
for tax in env_sign_sp:
    # Remove Dog others and Dog Shelter
    tax_w_metadata = tax_w_metadata[~tax_w_metadata['env_classification'].isin(['Dog Shelter', 'Dog others'])]
    categories = ['Dog Pet', 'Dog Colony', 'Dog Free_roaming', 'Wild Canid']
    custom_palette = sns.color_palette("Dark2", len(categories))

    # Figure size
    width_mm = 57
    height_mm = 35
    figsize_inch = (width_mm / 25.4, height_mm / 25.4)

    # Plot figure
    fig, ax = plt.subplots(figsize=figsize_inch)
    # Plot boxplot + swarmplot
    sns.swarmplot(data=tax_w_metadata, x='env_classification', y=tax,
                  order=categories, palette=custom_palette,
                  alpha=0.7, size=2, legend=False)
    sns.boxplot(data=tax_w_metadata, x='env_classification', y=tax,
                color='white', width=0.5, order=categories, showfliers=False)

    ax.set_title(tax,size=10, fontstyle='italic')
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x:.1f}'))
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_xlabel('')
    ax.set_ylabel('Rel ab.')

    sns.despine()
    plt.tight_layout()
    plt.show()
    out_path = 'intermediate-outputs/figures/maaslin2_tax/'+ tax + '_57env.svg'
    #plt.savefig(out_path)

### Age_classification: Plot boxplots and add pairwise comparisons
for tax in age_sign_sp:
    tax_w_metadata = tax_w_metadata[~tax_w_metadata['Animal_age_simplified'].isin(['unknown'])]
    categories = ['Young', 'Adult', 'Senior']
    custom_palette = {"Young": "#e6ab02", "Adult": "#d95f02", "Senior": "#a6761d"}

    # Figure size
    width_mm = 57
    height_mm = 35
    figsize_inch = (width_mm / 25.4, height_mm / 25.4)

    # Plot figure
    fig, ax = plt.subplots(figsize=figsize_inch)
    # Plot boxplot + swarmplot
    sns.swarmplot(data=tax_w_metadata, x='Animal_age_simplified', y=tax,
                  order=categories, palette=custom_palette,
                  alpha=0.7, size=2, legend=False)
    sns.boxplot(data=tax_w_metadata, x='Animal_age_simplified', y=tax,
                color='white', width=0.5, order=categories, showfliers=False)

    ax.set_title(tax,size=10, fontstyle='italic')
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x:.1f}'))
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_xlabel('')
    ax.set_ylabel('Rel ab.')

    sns.despine()
    plt.tight_layout()
    #plt.show()
    out_path = 'intermediate-outputs/figures/maaslin2_tax/'+ tax + '_57age.svg'
    plt.savefig(out_path)

### Compute Mann-Whitney U pairwise tests ###
pairwise_results = []
variable = 'Animal_age_simplified' #'env_classification'
sign_sp = age_sign_sp #env_sign_sp

# Iterate through significant variables
for tax in sign_sp:
    # Get unique groups for pairwise comparisons
    groups = sorted(tax_w_metadata[variable].unique())

    # Perform Mann-Whitney U test for each pair of groups
    for g1, g2 in combinations(groups, 2):
        group1 = tax_w_metadata[tax][tax_w_metadata[variable] == g1]
        group2 = tax_w_metadata[tax][tax_w_metadata[variable] == g2]

        # Ensure both groups have >1 sample
        if len(group1) > 1 and len(group2) > 1:
            stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
            pairwise_results.append((tax, g1, g2, stat, p))

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