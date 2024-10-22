import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

os.chdir('/data/Projects/ShanghaiDogs')
plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

# Input files
data = pd.read_csv('intermediate-outputs/external_datasets_mappings/reads_mapped_shd.tsv', sep='\t')
metadata = pd.read_csv('data/ShanghaiDogsMetadata/REP_canid_metadata.csv',sep=',')

# Calculate aligned fractions
data.eval('fraction = aligned/total', inplace=True)
data.eval('fraction_sp = aligned_sp/total', inplace=True)

## Plotting

# Fig_size
width_mm = 160
height_mm = 80
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)
data = data.melt(id_vars=['sample', 'group'])
data['variable'] = data['variable'].map({'fraction': 'SHD total MAGs', 'fraction_sp': 'SHD species-level catalog'})
data['value'] *= 100
data.rename(columns={
            'value': 'Reads mapped (%)',
            'variable': 'Reference',
            }, inplace=True)

# Group re-naming & ordering
data = pd.merge(data, metadata[['Sample_id','env_classification']],left_on='sample',right_on='Sample_id')
data['group'] = data['group']+' '+data['env_classification']
data['group'] = data['group'].str.replace('Berlin Dog Pet', 'Pet dogs\n(Berlin)')
data['group'] = data['group'].str.replace('Yarlagadda Dog Pet', 'Pet dogs\n(South Africa)')
data['group'] = data['group'].str.replace('Nestlé Dog Colony', 'Colony dogs\n(USA)')
data['group'] = data['group'].str.replace('Yarlagadda Dog Free_roaming', 'Non-pet dogs\n(India and Laos)')
data['group'] = data['group'].str.replace('Yarlagadda Dog Shelter', 'Non-pet dogs\n(India and Laos)')

group_order = [
    'Pet dogs\n(Berlin)',
    'Pet dogs\n(South Africa)',
    'Colony dogs\n(USA)',
    'Non-pet dogs\n(India and Laos)']

sns.boxplot(data=data, x='group', y='Reads mapped (%)', hue='Reference',
            order=group_order, ax=ax, boxprops=dict(alpha=.3), showfliers=False)
sns.stripplot(data=data, x='group', y='Reads mapped (%)', hue='Reference',
              order=group_order, dodge=True, ax=ax, legend=False)

ax.set_xlabel('')
ax.tick_params(labelsize=10)
ax.set_ylabel("Reads mapped (%)",fontsize=10)

# Capture the legend object
legend = ax.legend(loc='lower left', title="Reference", title_fontsize=10)
legend._legend_box.align = "left"

sns.despine(fig)
fig.tight_layout()

#plt.show()
fig.savefig('analysis/figures/reads_mapped_shd.svg')



