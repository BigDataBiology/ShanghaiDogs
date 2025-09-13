# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

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
dark2_custom_ord = ['#a6761d','#66a61e','#d95f02','#7570b3','#e7298a','#1b9e77','#666666']
phylum_to_color = {phylum: dark2_custom_ord[i % len(dark2_custom_ord)] for i, phylum in enumerate(phylum_counts_merged['Phylum'].unique())}
phylum_colors = phylum_counts_merged['Phylum'].map(phylum_to_color) # Assign colors based on Phylum

### Plotting donutplot for phyla distribution
width_mm = 65
height_mm = 90
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)
ax.set_position([0.05, 0.05, 0.8, 0.8])  # Adjust this for better centering
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
fig.savefig('intermediate-outputs/figures/donutplot_phyla_ver.svg')

### Plot stacked barplot, instead of donutplot

phylum_counts_merged['Percentage'] = (phylum_counts_merged['Count'] / phylum_counts_merged['Count'].sum()) * 100

width_mm = 60
height_mm = 130
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

# Create the figure and axes
fig, ax = plt.subplots(figsize=figsize_inch)

# Create the stacked bar
cumulative_height = 0
for i, (percentage, phylum, color) in enumerate(zip(
        phylum_counts_merged['Percentage'], phylum_counts_merged['Phylum'], phylum_colors
)):
    # Add each segment of the bar
    ax.bar(
        x=[0],
        height=[percentage],
        width=0.8,
        bottom=cumulative_height,
        color=[color],
        edgecolor='black'
    )

    # Calculate the middle of the segment for label placement
    middle_height = cumulative_height + percentage / 2

    # Add the label
    ax.text(
        x=0,
        y=middle_height,
        s=f"{phylum}\n({percentage:.1f}%)",
        ha='center',
        va='center',
        fontsize=9,
        color='black'
    )

    # Update cumulative height
    cumulative_height += percentage

# Adjust the plot
ax.set_xlim(-0.6, 0.6)
ax.set_ylim(0, cumulative_height)
ax.axis('off')  # Remove axes

# Tight layout and display
plt.tight_layout()
#plt.show()
plt.savefig('intermediate-outputs/figures/barplot_phyla.svg')

### Plot most prevalent species in SHD
sp_MAGs_counts = MIMAG_report['Species'].value_counts().reset_index()
sp_MAGs_counts_prev = sp_MAGs_counts[sp_MAGs_counts['count']>=30] # Genome assembled >30 times (~30 dogs)

# Link species to phylum info for ordering and coloring
sp_MAGs_counts_prev_phylum = (pd.merge(sp_MAGs_counts_prev,MIMAG_report[['Species','Phylum']],left_on='Species',right_on='Species')
                              .drop_duplicates(subset='Species'))
order = list(phylum_counts.index)

sp_MAGs_counts_prev_phylum['Phylum'] = pd.Categorical(
    sp_MAGs_counts_prev_phylum['Phylum'],
    categories=order,
    ordered=True
)

sp_MAGs_counts_prev_phylum = sp_MAGs_counts_prev_phylum.sort_values(
    by=['Phylum', 'Species'],
    ascending=[False, False]
)

# Calculate prevalence %
sp_MAGs_counts_prev_phylum['Percent']=sp_MAGs_counts_prev_phylum['count']/52*100
sp_MAGs_counts_prev_phylum = sp_MAGs_counts_prev_phylum.set_index('Species')

# Define a default color for phyla not in the dictionary
default_color = '#666666'  # Other phyla category
bar_colors = sp_MAGs_counts_prev_phylum['Phylum'].map(lambda phylum: phylum_to_color.get(phylum, default_color))

# Plot prevalence bar
width_mm = 90
height_mm = 120
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)
sp_MAGs_counts_prev_phylum['Percent'].plot(
    kind='barh',
    stacked=True,
    color=bar_colors,
    #alpha=0.8,
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
plt.savefig('intermediate-outputs/figures/high_prev_species.svg')

### Calculate abundant MAGs (for merging to previous plot)
median_RA_MAGs = repbin_cov.T.apply(lambda x: x[x != 0].median()).reset_index() # if present, calculate the median RA - ignore 0s
abundant_MAGs_sp = pd.merge(median_RA_MAGs,MIMAG_report['Species'],left_on='index',right_index=True)

prev_sp_ls = list(sp_MAGs_counts_prev_phylum.index)
prev_ab_MAGs_sp = abundant_MAGs_sp[abundant_MAGs_sp['Species'].isin(prev_sp_ls)]
prev_ab_MAGs_sp_phylum = pd.merge(prev_ab_MAGs_sp,MIMAG_report[['Species','Phylum']],right_on='Species',left_on='Species')
prev_ab_MAGs_sp_phylum = prev_ab_MAGs_sp_phylum.drop_duplicates('Species')

prev_ab_MAGs_sp_phylum['Phylum'] = pd.Categorical(
    prev_ab_MAGs_sp_phylum['Phylum'],
    categories=order,
    ordered=True)

prev_ab_MAGs_sp_phylum = prev_ab_MAGs_sp_phylum.sort_values(by=['Phylum', 'Species'],
                                                            ascending=[False, False])

prev_ab_MAGs_sp_phylum['Percent'] = prev_ab_MAGs_sp_phylum[0]*100
prev_ab_MAGs_sp_phylum = prev_ab_MAGs_sp_phylum.set_index('Species')

# Plot Bubbleplot
width_mm = 80
height_mm = 125
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)
scatter = ax.scatter(
    [1] * len(prev_ab_MAGs_sp_phylum),
    prev_ab_MAGs_sp_phylum.index,      # Index names on the y-axis
    s=prev_ab_MAGs_sp_phylum['Percent'] * 25,  # Scale bubble size
    alpha=0.9, edgecolors="#DDDDDD", c=bar_colors)

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
plt.savefig('intermediate-outputs/figures/high_prev-ab_sp_bubble_by_phyla.svg')

### Plot abundances distribution of most prevalent MAGs
ls_prev_ab_MAGs_sp = list(prev_ab_MAGs_sp_phylum.index)
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
CM.savefig('intermediate-outputs/figures/clustermap_ab_prevalent_species.svg')


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
