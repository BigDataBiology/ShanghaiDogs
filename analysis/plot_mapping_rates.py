import os
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

# Data prep for plotting
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
data['group'] = data['group'].str.replace('Berlin Dog Pet', 'Pet dogs (Berlin)')
data['group'] = data['group'].str.replace('Yarlagadda Dog Pet', 'Pet dogs (South Africa)')
data['group'] = data['group'].str.replace('Nestlé Dog Colony', 'Colony dogs (USA)')
data['group'] = data['group'].str.replace('Yarlagadda Dog Free_roaming', 'Non-pet dogs (India and Laos)')
data['group'] = data['group'].str.replace('Yarlagadda Dog Shelter', 'Non-pet dogs (India and Laos)')
data['group'] = data['group'].str.replace('NomNomNow Dog Pet', 'Pet dogs (USA)')
data['group'] = data['group'].str.replace('Shanghai Dog Pet', 'Pet dogs (Shanghai)')

group_order = [
    'Pet dogs (Shanghai)',
    'Pet dogs (Berlin)',
    'Pet dogs (South Africa)',
    'Pet dogs (USA)',
    'Colony dogs (USA)',
    'Non-pet dogs (India and Laos)']

### Plotting: VERTICAL boxplot

width_mm = 160
height_mm = 80
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)

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
fig.savefig('analysis/figures/reads_mapped_shd_vertical.svg')

### Plotting: HORIZONTAL boxplot

# Fig_size
width_mm = 155
height_mm = 50
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)

sns.boxplot(data=data, x='Reads mapped (%)', y='group', hue='Reference',
            order=group_order, ax=ax, boxprops=dict(alpha=.3), showfliers=False)
sns.stripplot(data=data, x='Reads mapped (%)', y='group', hue='Reference',
              order=group_order, size=2, dodge=True, ax=ax, legend=False)

ax.set_xlabel('Reads mapped (%)',fontsize=10)
ax.set_ylabel('')
ax.tick_params(labelsize=10)

# Add the legend at the figure level
ax.get_legend().remove()
fig.legend(
    loc='lower left',
    bbox_to_anchor=(0.01, 0.01),  # Adjust coordinates as needed
    title="")

sns.despine(fig)
fig.tight_layout()

#plt.show()
fig.savefig('analysis/figures/reads_mapped_shd_horizontal.svg')

### Plotting: HORIZONTAL boxplot w/o species-level MAGs

data_filt = data[data['Reference'] == 'SHD total MAGs']
data_filt = data_filt[data_filt['group'].isin(group_order)]

data_filt['env_classification']=data_filt['env_classification'].str.replace('Dog Free_roaming','Non-pet dog')
data_filt['env_classification']=data_filt['env_classification'].str.replace('Dog Shelter','Non-pet dog')

palette = ['#1b9e77', '#d95f02', '#7570b3'] #Dark2 colors

# Fig_size
plt.clf()

width_mm = 120
height_mm = 45
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)

sns.boxplot(data=data_filt, x='Reads mapped (%)', y='group', hue='env_classification',
            order=group_order, palette=palette, ax=ax, dodge=False, boxprops=dict(alpha=.3),
            showfliers=False, width=0.7)
sns.stripplot(data=data_filt, x='Reads mapped (%)', y='group', hue='env_classification',
              order=group_order, palette=palette, size=2, dodge=False, ax=ax, legend=False)

ax.set_xlabel('Reads mapped (%)',fontsize=10)
ax.set_ylabel('')
ax.tick_params(labelsize=10)
ax.get_legend().remove()

sns.despine(fig)
fig.tight_layout()

#plt.show()
fig.savefig('analysis/figures/reads_mapped_shd_by_env.svg')

