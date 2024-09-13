import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

data = pd.read_csv('../intermediate-outputs/external_datasets_mappings/reads_mapped_shd.tsv', sep='\t')
data.eval('fraction = aligned/total', inplace=True)
data.eval('fraction_sp = aligned_sp/total', inplace=True)
fig, ax = plt.subplots()
data = data.melt(id_vars=['sample', 'group'])
data['variable'] = data['variable'].map({'fraction': 'SHD', 'fraction_sp': 'SHD-Reps'})
data['value'] *= 100
data.rename(columns={
            'value': 'Reads mapped (%)',
            'variable': 'Reference',
            }, inplace=True)
sns.boxplot(data=data, x='group', y='Reads mapped (%)', hue='Reference', ax=ax, boxprops=dict(alpha=.3), showfliers=False)
sns.stripplot(data=data, x='group', y='Reads mapped (%)', hue='Reference', dodge=True, ax=ax, legend=False)
sns.despine(fig, trim=True)
fig.tight_layout()
fig.savefig('figures/reads_mapped_shd.pdf')



