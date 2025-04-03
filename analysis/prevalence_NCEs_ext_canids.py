import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.multitest import multipletests
from itertools import combinations

os.chdir('/data/Projects/ShanghaiDogs')
plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

# IMPORT INPUT FILES
metadata = pd.read_csv('data/ShanghaiDogsMetadata/REP_canid_metadata.csv', \
                       sep=',', encoding= 'unicode_escape', index_col=0)
metadata['Sex'] = metadata['Sex'].str.lower()

MAGs_NCE_covered_frac = pd.read_csv('intermediate-outputs/external_datasets_mappings/SHD_covered_fraction.tsv.gz',
                                    sep='\t', index_col=0, compression='gzip')
NCE_info = pd.read_csv('data/ShanghaiDogsTables/SHD_NC_props.tsv.gz',
                       sep='\t', index_col=0, compression='gzip')

# Prevalence estimates for NCEs
NCE_covered_frac = MAGs_NCE_covered_frac[MAGs_NCE_covered_frac.index.str.contains('NC.')]
NCE_covered_frac_prev = NCE_covered_frac.applymap(lambda x: 1 if x >= 0.8 else 0)
NCE_covered_frac_prev_T = NCE_covered_frac_prev.T
NCE_covered_frac_prev_T.index = NCE_covered_frac_prev_T.index.str.split('_').str[0]

# Table preparation for stats
NCE_covered_frac_metadata = pd.merge(NCE_covered_frac_prev_T,metadata[['Study','Size_class','Animal_age_simplified',
                                                                       'Sex','Body_condition','env_classification']],
                                     left_index=True,right_index=True)

metadata_cols = ['Study','Size_class','Animal_age_simplified',
                 'Sex','Body_condition','env_classification']

# Assess statistical significance
def analyze_NCE_assoc(df, metadata_cols, alpha=0.05):
    """
    Analyze associations between NCE presence/absence and metadata variables using Chi-Square
    and Fisher’s Exact tests for pairwise comparisons.
    Applies multiple testing correction across all NCEs and metadata variables (conservative, strict approach)

    Parameters:
    df (pd.DataFrame): DataFrame where rows contain samples, and columns contain NCEs + metadata variables.
    metadata_cols (list): List of column names that are metadata categories.
    alpha (float): Significance threshold for Chi-Square test (default = 0.05).

    Returns: a results DataFrame with globally corrected p-values.
    """

    results = []
    pairwise_pvals = []
    comparisons = []

    # Identify NCE columns (all except metadata columns)
    nce_cols = [col for col in df.columns if col not in metadata_cols]

    for nce in nce_cols:
        for metadata in metadata_cols:
            # Create contingency table
            contingency_table = pd.crosstab(df[nce], df[metadata])

            # Perform Chi-Square test
            chi2, p_chi2, _, _ = chi2_contingency(contingency_table)

            # Proceed to post-hoc tests only if Chi-Square is significant
            if p_chi2 < alpha:
                # Get unique categories in metadata variable
                categories = df[metadata].dropna().unique()

                # Pairwise Fisher’s Exact tests for each category pair
                for cat1, cat2 in combinations(categories, 2):
                    sub_df = df[df[metadata].isin([cat1, cat2])]
                    table = pd.crosstab(sub_df[nce], sub_df[metadata])

                    if table.shape == (2, 2):  # Fisher's Exact Test only for 2x2 tables
                        _, fisher_p = fisher_exact(table)
                        pairwise_pvals.append(fisher_p)
                        comparisons.append({"NCE": nce, "Comparison": f"{metadata}: {cat1} vs {cat2}"})

    # Apply Global Multiple Testing Correction (FDR across ALL NCEs and metadata)
    if pairwise_pvals:
        corrected_pvals = multipletests(pairwise_pvals, method='fdr_bh')[1]

        results = [
            {"NCE": comp["NCE"], "Comparison": comp["Comparison"], "Corrected_p": p_corr}
            for comp, p_corr in zip(comparisons, corrected_pvals)
        ]

    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    return results_df

# Run the analysis
results_df = analyze_NCE_assoc(NCE_covered_frac_metadata, metadata_cols)
results_df['Corrected_p_full'] = results_df['Corrected_p'].apply(lambda x: f"{x:.2e}")
results_df['Metadata'] = results_df['Comparison'].str.split(':').str[0]
results_df.to_csv('intermediate-outputs/tables/NCEs_comparisons_canid.csv',sep=',')

# Significant results
sign_results = results_df[results_df["Corrected_p"] < 0.05]
sign_results = sign_results[~sign_results["Comparison"].str.contains("unknown")]
sign_results = sign_results[~sign_results["Comparison"].str.contains("unclass")]
sign_nce_ls = list(sign_results['NCE'].unique())

# Create a dictionary with NCEs and the most significant metadata variable
sign_metadata_dict = {
    nce: sign_results.loc[sign_results['NCE'] == nce].nsmallest(1, 'Corrected_p')['Metadata'].tolist()
    for nce in sign_nce_ls}

# Filter prevalence table for significant NCE & representative samples
common_samples_ls = metadata.index.intersection(NCE_covered_frac_prev_T.index) # Representative samples + studies chosen
NCE_covered_frac_prev_filt = NCE_covered_frac_prev_T.loc[common_samples_ls, sign_nce_ls]

# Link to metadata
NCE_covered_frac_metadata = pd.merge(NCE_covered_frac_prev_filt,metadata,left_index=True,right_index=True)

# For each significant NCE assess their distribution according significant metadata variables
print(sign_metadata_dict)

# Create stacked bar plots
width_mm = 45
height_mm = 40
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

category_orders = {
    'env_classification': ['Dog Pet', 'Dog Colony', 'Dog Shelter', 'Dog Free_roaming'],
    'Size': ['small', 'medium', 'large'],
    'Sex': ['male', 'female'],
    'Animal_age_simplified': ['Young', 'Adult', 'Senior'],
    'Study': ['This_study', 'Berlin_cohort', 'Tanprasertsuk_2021', 'Coelho_2018',
              'Allaway_2020_dogs','Yarlagadda_2022_global_dog'],
    'Body_condition': ['Lean to normal','Overweight to obese']}

abbreviations = {
    'Dog Pet': 'P',
    'Dog Colony': 'C',
    'Dog Shelter': 'S',
    'Dog Free_roaming': 'F',
    'Small': 'S',
    'Medium': 'M',
    'Large': 'L',
    'male': 'M',
    'female': 'F',
    'Young': 'Y',
    'Adult': 'A',
    'Senior': 'S',
    'Lean to normal': 'N',
    'Overweight to obese': 'O',
    'This_study': 'SHD',
    'Berlin_cohort': 'BER',
    'Tanprasertsuk_2021': 'US_P',
    'Coelho_2018': 'US_C',
    'Allaway_2020_dogs':'UK_C',
    'Yarlagadda_2022_global_dog':'YAR'
}

#nce = 'SHD1_NC.011'
#meta = 'env_classification'

for nce in sign_nce_ls:
    for meta in sign_metadata_dict[nce]:
        plt.figure(figsize=figsize_inch)

        if meta in category_orders:
            NCE_covered_frac_metadata[meta] = pd.Categorical(
                NCE_covered_frac_metadata[meta],
                categories=category_orders[meta],
                ordered=True)

        # Count occurrences of NCE presence by metadata category
        grouped_data = NCE_covered_frac_metadata.groupby(meta)[nce].value_counts(normalize=True).unstack().fillna(0)
        grouped_data = grouped_data[[1, 0]]

        # Create stacked bar plot
        ax = grouped_data.plot(kind='bar', stacked=True, color=['#d95f02','#aaaaaa'],
                               width=0.8, ax=plt.gca())

        # Update x-axis labels using abbreviations
        new_labels = [abbreviations.get(label, label) for label in grouped_data.index]  # Replace if abbreviation exists
        ax.set_xticklabels(new_labels, ha='center',rotation=0)

        # Formatting
        sns.despine()
        plt.title(f'{nce}')
        plt.ylabel('')
        plt.xlabel('')
        plt.legend().remove()
        #plt.xticks([], [])
        plt.tight_layout()

        # Show plot
        # plt.show()
        out_path = 'intermediate-outputs/figures/non-chr-elements/' + nce + '_' + meta + '.svg'
        print('Saving figure at ' + out_path)
        plt.savefig(out_path)


# Plot clustermap - too crowded
width_mm = 800
height_mm = 400
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

# Set up the clustermap directly (since it doesn’t use fig, ax format natively)
CM = sns.clustermap(NCE_covered_frac_prev_filt, figsize=figsize_inch, row_cluster=True,
                    cmap='Blues', linewidth=0.3,
                    yticklabels=False, xticklabels=True)
CM.ax_heatmap.set_ylabel('')
CM.ax_heatmap.tick_params(axis='x', labelsize=10)

plt.show()
