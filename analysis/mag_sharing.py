import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

data = pd.read_csv('../intermediate-outputs/05_dereplication/01_drep/ANI_9999/data_tables/Ndb.csv')
metadata = pd.read_csv('../data/ShanghaiDogsMetadata/SH_Dog_metadata_red.csv', index_col=0)
household = metadata['Household.1'].to_dict()
household['D000'] = household['D008']

data['ref_sample'] = data.reference.str.split('.').str[0].str.split('_').str[-1]
data['qry_sample'] = data.querry.str.split('.').str[0].str.split('_').str[-1]
data = data[data['alignment_coverage']>0.50]

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
fig.savefig('figures/mag_sharing.pdf')

print(f'''Mann-Whitney U test: p-value = {m_v.pvalue}''')

