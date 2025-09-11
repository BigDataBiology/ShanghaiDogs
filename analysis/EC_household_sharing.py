import os
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations
from sklearn.metrics import jaccard_score
import seaborn as sns
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

os.chdir('/data/Projects/ShanghaiDogs')
plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

# Input data
mimag_tab = pd.read_csv('data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv',index_col=0)
metadata = pd.read_csv('data/ShanghaiDogsMetadata/SH_Dog_metadata_red.csv', index_col=0)
metadata = metadata[~metadata.index.isin(['D000'])] # remove biological replicate

MAGs_NCE_covered_frac = pd.read_csv('intermediate-outputs/external_datasets_mappings/SHD_covered_fraction.tsv.gz',
                                    sep='\t', index_col=0, compression='gzip')
NCE_info = pd.read_csv('data/ShanghaiDogsTables/SHD1_EC_props.tsv.gz',
                       sep='\t', index_col=0, compression='gzip')

# Prevalence estimates for NCEs
NCE_covered_frac = MAGs_NCE_covered_frac[MAGs_NCE_covered_frac.index.str.contains('NC.')]
NCE_covered_frac_prev = NCE_covered_frac.applymap(lambda x: 1 if x >= 0.8 else 0)
NCE_covered_frac_prev_T = NCE_covered_frac_prev.T
NCE_covered_frac_prev_T.index = NCE_covered_frac_prev_T.index.str.split('_').str[0]

# Keep only SHD samples
NCE_covered_frac_prev_T = NCE_covered_frac_prev_T[NCE_covered_frac_prev_T.index.str.contains('D0')]
NCE_covered_frac_prev_T = NCE_covered_frac_prev_T[~NCE_covered_frac_prev_T.index.isin(['D000'])]

# Split prevalence table by NCE type
NCE_plasmid_ls = NCE_info[NCE_info['Category']=='plasmid'].index.to_list()
NCE_virus_ls = NCE_info[NCE_info['Category']=='virus'].index.to_list()
NCE_uncat_ls = NCE_info[NCE_info['Category']=='uncategorized'].index.to_list()

NCE_plasmid = NCE_covered_frac_prev_T.loc[:, NCE_covered_frac_prev_T.columns.isin(NCE_plasmid_ls)]
NCE_virus = NCE_covered_frac_prev_T.loc[:, NCE_covered_frac_prev_T.columns.isin(NCE_virus_ls)]
NCE_uncat = NCE_covered_frac_prev_T.loc[:, NCE_covered_frac_prev_T.columns.isin(NCE_uncat_ls)]

# Filter for minimum sample prevalence
def filter_by_prevalence(df, min_prevalence=0.1):
    """
    Keep only columns (ECs) that are present in at least `min_prevalence` of samples (index).
    min_prevalence is a proportion
    """
    prevalence = df.sum(axis=0) / df.shape[0]
    return df.loc[:, prevalence >= min_prevalence]

NCE_plasmid_filt = filter_by_prevalence(NCE_plasmid)
NCE_virus_filt = filter_by_prevalence(NCE_virus)
NCE_uncat_filt = filter_by_prevalence(NCE_uncat)
NCE_all_filt = filter_by_prevalence(NCE_covered_frac_prev_T)

# Compute Jaccard similarity
def compute_pairwise_jaccard(df, sample_metadata, ec_type):
    results = []
    for s1, s2 in combinations(df.index, 2):
        sim = jaccard_score(df.loc[s1], df.loc[s2],zero_division=0)
        same_house = sample_metadata.loc[s1, 'Household.1'] == sample_metadata.loc[s2, 'Household.1']
        results.append({
            'sample1': s1,
            'sample2': s2,
            'jaccard_similarity': sim,
            'same_household': same_house,
            'ec_type': ec_type
        })
    return pd.DataFrame(results)

Jacc_all = compute_pairwise_jaccard(NCE_all_filt,metadata,'all')
Jacc_plasmid = compute_pairwise_jaccard(NCE_plasmid_filt,metadata,'plasmid')
Jacc_virus = compute_pairwise_jaccard(NCE_virus_filt,metadata,'virus')
Jacc_uncat = compute_pairwise_jaccard(NCE_uncat_filt,metadata,'uncategorized')

# Concatenate the NCEs results
jacc_all_by_type = pd.concat([Jacc_plasmid,Jacc_virus,Jacc_uncat], ignore_index=True)

# Compute statistics
stat_results = []
jacc_df = Jacc_all #Jacc_all jacc_all_by_type

for ec in jacc_df['ec_type'].unique():
    data = jacc_df[jacc_df['ec_type'] == ec]
    same = data[data['same_household'] == True]['jaccard_similarity']
    diff = data[data['same_household'] == False]['jaccard_similarity']

    # Mann-Whitney U test
    stat, pval = mannwhitneyu(same, diff, alternative='greater')  # test if same > diff

    stat_results.append({
        'ec_type': ec,
        'n_same_house': len(same),
        'n_diff_house': len(diff),
        'mean_same_house': same.mean(),
        'mean_diff_house': diff.mean(),
        'statistic': stat,
        'p_value': pval
    })

results_df = pd.DataFrame(stat_results)
results_df['p_adj'] = multipletests(results_df['p_value'], method='fdr_bh')[1]

# Visualization jacc_all_by_type
fig, ax = plt.subplots(figsize=(3, 3))  # adjust size as needed

sns.boxplot(data=jacc_all_by_type, x='ec_type', y='jaccard_similarity',
            width=0.7, ax=ax, palette='Dark2', hue='same_household',
            flierprops={'markersize': 3})
ax.set_ylabel('Jaccard Similarity', fontsize=10)
ax.set_xlabel('', fontsize=10)
ax.set_title('Extra-chromosomal elements', fontsize=11)
ax.set_xticklabels(['plasmid', 'virus','uncat'],rotation=45,ha='right')
ax.legend_.remove()

# Tight layout and show
sns.despine()
plt.tight_layout()
plt.show()

# Visualization jacc_all
width_mm = 75
height_mm = 75
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)  # adjust size as needed

#sns.boxplot(data=Jacc_all, x='same_household', y='jaccard_similarity',
#            width=0.7, ax=ax, palette='Dark2', flierprops={'markersize': 3})
ax.clear()
sns.boxplot(data=Jacc_all, x='same_household', y='jaccard_similarity', ax=ax,
            boxprops={'facecolor':'None'}, showfliers=False, width=0.7)
sns.stripplot(data=Jacc_all, x='same_household', y='jaccard_similarity', ax=ax,
              alpha=0.5, size=2.5, jitter=0.2, palette=['#D3D3D3','#e7298a'])

ax.set_title('Extra-chromosomal elements', fontsize=11, pad=22)
ax.text(0.5, 1.1, f'p-value = {pval:.1e}', ha='center', va='center', transform=ax.transAxes)
ax.set_ylim(0,1)

ax.set_ylabel('Jaccard Similarity', fontsize=10)
ax.set_xlabel('', fontsize=10)
ax.set_xticklabels(['Different\nHousehold', 'Same\nHousehold'],rotation=0,ha='center')

# Tight layout and show
sns.despine()
plt.tight_layout()
#plt.show()
fig.savefig('intermediate-outputs/figures/EC_sharing_Jacc.svg')
