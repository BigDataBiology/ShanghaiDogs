import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import numpy as np

os.chdir('/data/Projects/ShanghaiDogs/analysis')
data = pd.read_csv('../intermediate-outputs/05_dereplication/01_drep/ANI_9999/data_tables/Ndb.csv')
mimag_tab = pd.read_csv('../data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv',index_col=0)
metadata = pd.read_csv('../data/ShanghaiDogsMetadata/SH_Dog_metadata_red.csv', index_col=0)
household = metadata['Household.1'].to_dict()
household['D000'] = household['D008']

data['ref_sample'] = data.reference.str.split('.').str[0].str.split('_').str[-1]
data['qry_sample'] = data.querry.str.split('.').str[0].str.split('_').str[-1]
data = data[data['alignment_coverage']>0.50]
data = data[data['ani']>0.95]

paired = data.groupby(['ref_sample', 'qry_sample'])['ani'].apply(list)
paired = paired.reset_index()
paired['shared_sp'] = paired['ani'].apply(lambda x: sum(x > 0.95 for x in x))
paired['shared_st'] = paired['ani'].apply(lambda x: sum(x > 0.99 for x in x))
paired = paired[['ref_sample', 'qry_sample', 'shared_sp', 'shared_st']]
paired.query('ref_sample == "D000" and qry_sample == "D008"')
paired.eval('fraction_shared_st = shared_st/shared_sp', inplace=True)
paired['is_same_household'] = paired[['ref_sample', 'qry_sample']].apply(lambda x: household[x[0]] == household[x[1]], axis=1)
paired = paired.query('ref_sample < qry_sample')

paired_10 = paired.query('shared_sp >= 10')
m_v = stats.mannwhitneyu(paired_10.query('is_same_household')['fraction_shared_st'], paired_10.query('~is_same_household')['fraction_shared_st'])

fig, ax = plt.subplots()
ax.clear()
sns.boxplot(data=paired_10, x='is_same_household', y='fraction_shared_st', ax=ax, boxprops={'facecolor':'None'}, showfliers=False)
sns.stripplot(data=paired_10, x='is_same_household', y='fraction_shared_st', ax=ax)
ax.text(0.5, 0.99, f'Mann-Whitney U test: p-value = {m_v.pvalue:.1e}', ha='center', va='center', transform=ax.transAxes)
sns.despine(fig, trim=True)
fig.savefig('figures/mag_sharing_ANI.pdf')

print(f'''Mann-Whitney U test: p-value = {m_v.pvalue}''')

## Find bacterial species that are commonly >99% identical across SH dog cohort
## Find taxonomy of highly-shared MAGs

data['reference'] = data['reference'].str.replace('.fa.gz','')
data['querry'] = data['querry'].str.replace('.fa.gz','')
data_tax = pd.merge(data,mimag_tab[['Classification','Original ID']],left_on='reference',right_on='Original ID')
data_tax = pd.merge(data_tax,mimag_tab[['Classification','Original ID']],left_on='querry',right_on='Original ID')
data_tax.rename(columns={'Classification_x': 'Ref Classification', 'Classification_y': 'Qry Classification'}, inplace=True)
data_tax.drop(['Original ID_x','Original ID_y'],axis=1,inplace=True)

# Remove comparisons with eachself & for comparisons in the two directions, keep one
# Remove D000 comparisons, since D008 is the same
mask = data_tax['reference'] != data_tax['querry']
data_tax = data_tax[mask]
data_tax['Comparison'] = data_tax[['reference', 'querry']].apply(lambda x: tuple(sorted(x)), axis=1)
data_tax = data_tax.drop_duplicates('Comparison')
data_tax.drop(columns=['Comparison'], inplace=True)
data_tax = data_tax[~(data_tax['reference'].str.contains('D000'))]
data_tax = data_tax[~(data_tax['querry'].str.contains('D000'))]

# Keep within species comparisons only (according to GTDB-tk)
diffs = data_tax['Ref Classification'] != data_tax['Qry Classification']
diffs = data_tax[diffs]
print(diffs)

same_sp = data_tax['Ref Classification'] == data_tax['Qry Classification']
data_tax_sp = data_tax[same_sp]

# Proportion of same strain within the same species
data_tax_99 = data_tax_sp[data_tax_sp['ani']>0.99]
tax_grouping_99 = data_tax_99['Ref Classification'].value_counts().reset_index()
tax_grouping = data_tax['Ref Classification'].value_counts().reset_index()

merged_tax = pd.merge(tax_grouping_99,tax_grouping,right_on='Ref Classification',left_on='Ref Classification')
merged_tax = merged_tax[merged_tax['count_y'] > 100]
merged_tax['comparisons >99 ANI'] = merged_tax['count_x']/merged_tax['count_y']*100
merged_tax.columns = ['Classification','ani >=99%','ani >=95%','% ani >=99% if same sp']

merged_tax[['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']] = \
    merged_tax['Classification'].str.split(';', expand=True)

merged_tax_sorted = merged_tax.sort_values(by='% ani >=99% if same sp', ascending=False)
high_st_sharing = merged_tax_sorted[merged_tax_sorted['% ani >=99% if same sp']>70]
low_st_sharing = merged_tax_sorted[merged_tax_sorted['% ani >=99% if same sp']<10]

plt.figure(figsize=(8, 4))
order_colors = {phylum: plt.cm.Dark2(i) for i, phylum in enumerate(merged_tax_sorted['Phylum'].unique())}
plt.bar(high_st_sharing['Species'], high_st_sharing['% ani >=99% if same sp'], \
        color=[order_colors[order] for order in high_st_sharing['Phylum']])
plt.xticks(rotation=30, ha='right')
plt.ylim(0,100)
plt.tight_layout()
#plt.show()
plt.savefig('figures/species_high_ANI.pdf')
