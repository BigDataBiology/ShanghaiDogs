# -*- coding: utf-8 -*-

import os
import pandas as pd
from scipy.stats import trim_mean
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.cm import get_cmap
from matplotlib.colors import to_hex

plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot
os.chdir('/data/Projects/ShanghaiDogs')

### Import and format data
repbin_cov = pd.read_csv('intermediate-outputs/repbin_coverage_rmean.tsv',
                         index_col=0, sep='\t')
repbin_cov.columns = repbin_cov.columns.str.replace('_SR_to_95_ANI', '', regex=False) # remove "_SR_to_95_ANI" in column headers

MIMAG_report = pd.read_csv('data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv', \
                    delimiter=',', header=0, index_col=0)
MIMAG_report.index = MIMAG_report.index.str.replace('.fna.gz', '', regex=False) # remove ".fna.gz" in index names
MIMAG_report['Phylum'] = MIMAG_report['Classification'].str.extract(r'p__([^;]+)')
MIMAG_report['Family'] = MIMAG_report['Classification'].str.extract(r'f__([^;]+)')
MIMAG_report['Genus'] = MIMAG_report['Classification'].str.extract(r'g__([^;]+)')
MIMAG_report['Species'] = MIMAG_report['Classification'].str.extract(r's__([^;]+)')

### Total MAG count by phylum
phylum_counts = MIMAG_report['Phylum'].value_counts()
phylum_large = phylum_counts[phylum_counts >= 100] # At least 100 total MAGs from that phylum
phylum_small = phylum_counts[phylum_counts < 100]
phylum_counts_merged = phylum_large.copy()
phylum_counts_merged['Other phyla'] = phylum_small.sum()
phylum_counts_merged = phylum_counts_merged.reset_index()
phylum_counts_merged.columns = ['Phylum','Count']

# Color palette
dark2_colors = get_cmap("Dark2_r").colors
custom_palette = [to_hex(color) for i, color in enumerate(dark2_colors) if i not in [2]]  # Exclude green and yellow
phylum_to_color = {phylum: custom_palette[i % len(custom_palette)] for i, phylum in enumerate(phylum_counts_merged['Phylum'].unique())}
phylum_colors = phylum_counts_merged['Phylum'].map(phylum_to_color) # Assign colors based on Phylum

# Plotting donutplot
width_mm = 65
height_mm = 90
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

#colors = sns.color_palette("Dark2", len(phylum_counts_merged))

fig, ax = plt.subplots(figsize=figsize_inch)
ax.set_position([0.05, 0.05, 0.8, 0.8])  # Adjust this as needed for better centering
wedges, texts, autotexts = ax.pie(
    phylum_counts_merged['Count'],
    startangle=210,
    colors=phylum_colors,
    autopct='%1.1f%%',
    pctdistance=1.22,
    wedgeprops=dict(width=0.4, edgecolor='white')) # for donutplot

for autotext in autotexts:
    autotext.set_color('black')
    autotext.set_fontsize(10)

# Add legend
ax.legend(wedges, phylum_counts_merged['Phylum'], loc="upper center",
          bbox_to_anchor=(0.5,0.05), fontsize=10, ncol=1)

plt.tight_layout()
#plt.show()
fig.savefig('analysis/figures/donutplot_phyla_ver.svg')

# List of prevalent MAGs (>30 MAGs in SHD)
sp_MAGs_counts = MIMAG_report['Species'].value_counts().reset_index()
sp_MAGs_counts_prev = sp_MAGs_counts[sp_MAGs_counts['count']>=30] # Genome assembled >30 times (~30 dogs)

# order them according to MAG qual graph
order = ['Blautia sp000432195', 'Blautia_A caecimuris',
         'Fusobacterium_B sp900541465', 'Blautia_A sp900541345',
         'Blautia hansenii', 'Oliverpabstia sp000432335',
         'Collinsella intestinalis', 'Faecalimonas sp900550235',
         'Ruminococcus_B gnavus', 'Megamonas funiformis',
         'Schaedlerella glycyrrhizinilytica_A', 'Phocaeicola sp900546645',
         'Faecalibacterium sp900540455', 'Amedibacterium intestinale',
         'Enterocloster sp001517625', 'Peptacetobacter hiranonis',
         'Thomasclavelia spiroformis_A', 'Ventrimonas sp900538475',
         'Phocaeicola coprocola',  'Sutterella wadsworthensis_A',
         'Faecalimonas umbilicata', 'Clostridium_Q sp000435655',
         'Bacteroides sp900766005', 'Amedibacillus dolichus',
         'Eisenbergiella sp900539715']

sp_MAGs_counts_prev['Species'] = pd.Categorical(
    sp_MAGs_counts_prev['Species'],
    categories=order,
    ordered=True
)
sp_MAGs_counts_prev = sp_MAGs_counts_prev.sort_values(by='Species', ascending=True)
sp_MAGs_counts_prev['Percent']=sp_MAGs_counts_prev['count']/52*100
sp_MAGs_counts_prev = sp_MAGs_counts_prev.set_index('Species')

# Link species to phylum for coloring
MIMAG_report_unique_sp = MIMAG_report.drop_duplicates(subset='Species')
MIMAG_report_unique_sp = MIMAG_report_unique_sp.set_index('Species')
sp_MAG_prev = pd.merge(sp_MAGs_counts_prev,MIMAG_report_unique_sp['Phylum'],left_index=True,right_index=True)

# Define a default color for phyla not in the dictionary
default_color = '#1b9e77'  # Other phyla category
bar_colors = sp_MAG_prev['Phylum'].map(lambda phylum: phylum_to_color.get(phylum, default_color))

# Plot prevalence bar
width_mm = 90
height_mm = 110
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)
sp_MAG_prev['Percent'].plot(
    kind='barh',
    stacked=True,
    color=bar_colors,
    alpha=0.8,
    ax=ax,
    width=0.7
)

# Add a line at 50%
ax.axvline(50, color='black', linestyle='--', linewidth=0.7)

# Axis aesthetics
ax.set_xlabel('',fontsize=10)
ax.set_xticks([50, 100])
ax.set_xticklabels(['50%','100%'])
ax.set_ylabel('')
ax.set_xlim(35, 100)

plt.tight_layout()
#plt.show()
plt.savefig('analysis/figures/high_prev_species.svg')

# Calculate abundant MAGs
median_RA_MAGs = repbin_cov.T.apply(lambda x: x[x != 0].median()).reset_index() # if present, calculate the median RA - ignore 0s
abundant_MAGs_sp = pd.merge(median_RA_MAGs,MIMAG_report['Species'],left_on='index',right_index=True)
prev_ab_MAGs_sp = abundant_MAGs_sp[abundant_MAGs_sp['Species'].isin(order)]
prev_ab_MAGs_sp['Species'] = pd.Categorical(
    prev_ab_MAGs_sp['Species'],
    categories=order,
    ordered=True)
prev_ab_MAGs_sp = prev_ab_MAGs_sp.sort_values(by='Species', ascending=True)
prev_ab_MAGs_sp['Percent'] = prev_ab_MAGs_sp[0]*100
prev_ab_MAGs_sp = prev_ab_MAGs_sp.set_index('Species')

# Plot Bubbleplot
width_mm = 80
height_mm = 115
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)
scatter = ax.scatter(
    [1] * len(prev_ab_MAGs_sp),
    prev_ab_MAGs_sp.index,      # Index names on the y-axis
    s=prev_ab_MAGs_sp['Percent'] * 25,  # Scale bubble size
    alpha=0.8, edgecolors="#DDDDDD", c=bar_colors)

# Add legend for bubble sizes
bubble_sizes = [1, 5, 10]
for size in bubble_sizes:
    ax.scatter([], [], s=size * 25, color="lightgrey", edgecolors="#DDDDDD", label=f"{size}%")

fig.legend(loc="upper left", fontsize=10, title_fontsize=10)

# Axis aesthetics
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.tick_params(left = False, bottom = False)

# Show the plot
plt.tight_layout()
#plt.show()
plt.savefig('analysis/figures/high_prev-ab_sp_bubble_by_phyla__.svg')


# Plot abundances distribution of most prevalent MAGs
ls_prev_ab_MAGs_sp = list(prev_ab_MAGs_sp.index)
prev_ab_species = pd.merge(repbin_cov,MIMAG_report['Species'],left_index=True,right_index=True)
prev_ab_species = prev_ab_species.set_index('Species')
prev_ab_species = prev_ab_species[prev_ab_species.index.isin(ls_prev_ab_MAGs_sp)]

# Plot clustermap
width_mm = 80
height_mm = 170
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

# Convert data to percentages
prev_ab_species_percent = prev_ab_species * 100
prev_ab_species_percent = prev_ab_species_percent.T

# Set up the clustermap directly (since it doesn’t use fig, ax format natively)
CM = sns.clustermap(prev_ab_species_percent, figsize=figsize_inch, row_cluster=True,
                    cmap='YlOrBr', linewidth=0.3, vmax=5,
                    yticklabels=False, xticklabels=True)
CM.ax_heatmap.set_ylabel('')
CM.ax_heatmap.tick_params(axis='x', labelsize=10)

# Modify color bar to reflect percentage
cbar = CM.ax_heatmap.collections[0].colorbar
cbar.set_ticks([0, 1, 2, 3, 4, 5])
cbar.set_ticklabels([f"{int(t)}%" for t in cbar.get_ticks()])

# Adjust color bar position if needed
CM.ax_cbar.set_position((0.02, 0.78, 0.03, 0.2))

# Save or show plot
# plt.show()
CM.savefig('analysis/figures/clustermap_ab_prevalent_species.svg')


### Plot Heatmap
fig, ax = plt.subplots(figsize=figsize_inch)
sns.heatmap(prev_ab_species_percent, ax=ax, cmap='YlOrBr', linewidths=0.3,
            vmax=5, yticklabels=False, xticklabels=False, cbar_kws={'label': 'Percentage (%)'})
# Remove x and y labels
ax.set_xlabel('')
ax.set_ylabel('')

# Adjust color bar ticks to show percentages
cbar = ax.collections[0].colorbar
cbar.set_ticks([0, 1, 2, 3, 4, 5])
cbar.set_ticklabels([f"{int(t)}%" for t in cbar.get_ticks()])

plt.tight_layout()
plt.show()
