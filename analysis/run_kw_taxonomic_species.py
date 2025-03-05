"""
Created on Mon Mar 03 17:18:14 2025
@author: Anna Cusco
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests

os.chdir('/data/Projects/ShanghaiDogs/')
plt.rcParams['svg.fonttype'] = 'none'  # to avoid transforming the font to plot

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
tax_w_metadata = pd.merge(tax_profile_filt,
                          metadata_filt[['env_classification','Animal_age_simplified','Sex','Study','Size_class']],
                          left_index=True,right_index=True)

### Check KW diff between groups, as an alternative approach to Maaslin2
variable_ls = ['env_classification','Animal_age_simplified','Sex','Study','Size_class']
kruskal_all = []

# Filter data function
def filter_data(df, variable):
    study_ls = ['This_study', 'Coelho_2018', 'Allaway_2020_dogs', 'Bai_2023', 'Wang_2019_dogs']
    filtered = df.loc[~df[variable].isin(['Dog others', 'Dog Shelter', 'unclass', 'unknown', 'Unknown'])]
    if variable in ['Animal_age_simplified', 'Sex', 'Size_class']:
        filtered = filtered[filtered['env_classification'] != 'Wild Canid']
    if variable == 'Study':
        filtered = filtered[filtered[variable].isin(study_ls)]
    return filtered

# Loop through each variable and test across all taxonomies
for variable in variable_ls:
    print(f"Testing variable: {variable}")
    tax_w_metadata = pd.merge(tax_profile_filt, metadata_filt[
        ['env_classification', 'Animal_age_simplified', 'Sex', 'Study', 'Size_class']],
                              left_index=True, right_index=True)
    # Filter the data for the variable
    tax_metadata_filt = filter_data(tax_w_metadata, variable)
    print(tax_metadata_filt[variable].value_counts())

    # Loop through each bacterial species/taxon in the dataset
    for tax in tax_metadata_filt.columns:
        # Skip the non-taxonomic columns like metadata
        if tax in tax_metadata_filt.select_dtypes(include=['object', 'category']).columns:
            continue

        # Group the data by the categorical variable
        grouped_data = [tax_metadata_filt[tax_metadata_filt[variable] == g][tax].values
                        for g in tax_metadata_filt[variable].unique() if
                        len(tax_metadata_filt[tax_metadata_filt[variable] == g][tax].values) > 1]

        # Check for zero variance and exclude that group
        if len(grouped_data) > 1:
            group_variances = [np.var(group) for group in grouped_data]
            non_zero_variance_groups = [group for group, var in zip(grouped_data, group_variances) if var > 0]

            if len(non_zero_variance_groups) > 1:  # Proceed only if at least two groups have non-zero variance
                # Perform the Kruskal-Wallis test
                kw_stat, kw_pval = kruskal(*non_zero_variance_groups)

                # Store the results (variable, taxon, statistic, p-value)
                kruskal_all.append((variable, tax, kw_stat, kw_pval))

# Convert the results to a DataFrame for easy interpretation
kruskal_df = pd.DataFrame(kruskal_all, columns=['Variable', 'Taxonomy', 'KW_Statistic', 'p_value'])

# Bonferroni correction
kruskal_df['p_value_bonferroni'] = multipletests(kruskal_df['p_value'], method='bonferroni')[1]

# Benjamini-Hochberg FDR correction
kruskal_df['p_value_fdr'] = multipletests(kruskal_df['p_value'], method='fdr_bh')[1]

# Print or save the results
print(kruskal_df.head())

# Group by 'Taxonomy' and find the row with the lowest 'p_value_fdr' and highest 'KW_Statistic'
most_significant_kw = kruskal_df.loc[kruskal_df.groupby('Taxonomy')['p_value_fdr'].idxmin()]
most_significant_kw = most_significant_kw[most_significant_kw['p_value_fdr']<0.05]
most_significant_kw = most_significant_kw.sort_values(by='p_value_fdr', ascending=True)
most_significant_kw['p_value_fdr'] = most_significant_kw['p_value_fdr'].apply(lambda x: f"{x:.4e}")

env_sign = most_significant_kw[most_significant_kw['Variable']=='env_classification']
age_sign = most_significant_kw[most_significant_kw['Variable']=='Animal_age_simplified']
size_sign = most_significant_kw[most_significant_kw['Variable']=='Size_class']
study_sign = most_significant_kw[most_significant_kw['Variable']=='Study']
