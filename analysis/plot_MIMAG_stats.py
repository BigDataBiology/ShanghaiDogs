# -*- coding: utf-8 -*-

"""
Created on Mon Sep 02 10:41:58 2024
@author: Anna Cusco
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

### Import data
MIMAG_report = pd.read_csv('/data/Projects/ShanghaiDogs/data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv', \
                    delimiter=',', header=0, index_col=0)

MIMAG_report['Quality_det']=MIMAG_report['Quality']
for n in MIMAG_report.index:
    if MIMAG_report.loc[n, 'Nr contigs'] == 1 and MIMAG_report.loc[n, 'Quality'] == 'high-quality':
        MIMAG_report.loc[n, 'Quality_det'] = 'single-contig HQ'


### HQ_mq MAGs table
quality_table = MIMAG_report.pivot_table(
    index='Sample',
    columns='Quality',
    aggfunc='size',
    fill_value=0
)

quality_table['Total']=quality_table['high-quality']+quality_table['medium-quality']
quality_table = quality_table.sort_values('Total',ascending=False)

# Plot a stacked barplot
# Fig_size
width_mm = 110
height_mm = 60
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

# Barplot
ax = quality_table[['high-quality', 'medium-quality']].plot(kind="bar", stacked=True,color=['#1B9E77','#EDB458'],width=0.75,figsize=figsize_inch)

ax.set_ylabel("N# of MAGs",fontsize=10)
ax.set_xlabel('Samples',fontsize=10)
ax.set_yticklabels(map(int, ax.get_yticks()), fontsize=10)
ax.set_xticklabels([])
ax.tick_params(axis='x', bottom=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Custom legend text
legend_labels = ['High-quality MAGs', 'Medium-quality MAGs']
ax.legend(labels=legend_labels, title='', loc='upper right',fontsize=10)

# Save figure
plt.tight_layout()
# plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/analysis/figures/total_MAGs_qual.svg")


### HQ_mq 1-contig
quality_table_det = MIMAG_report.pivot_table(
    index='Sample',
    columns='Quality_det',
    aggfunc='size',
    fill_value=0
)

quality_table_det['Total']=quality_table_det['high-quality']+quality_table_det['medium-quality']+quality_table_det['single-contig HQ']
quality_table_det = quality_table_det.sort_values('Total',ascending=False)

# Plot a stacked barplot
# Fig_size
width_mm = 110
height_mm = 60
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

# Barplot
ax = (quality_table_det[['single-contig HQ','high-quality', 'medium-quality']].plot
      (kind='bar', stacked=True,color=['#1E3F20','#1B9E77','#EDB458'],width=0.75,
       figsize=figsize_inch))

ax.set_ylabel('N# of MAGs',fontsize=10)
ax.set_xlabel('Samples',fontsize=10)
ax.set_yticklabels(map(int, ax.get_yticks()), fontsize=10)
ax.set_xticklabels([])
ax.tick_params(axis='x', bottom=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Custom legend text
legend_labels = ['High-quality MAGs (1 contig)', 'High-quality MAGs', 'Medium-quality MAGs']
ax.legend(labels=legend_labels, title='', loc='upper right',fontsize=10)

# Save figure
plt.tight_layout()
# plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/analysis/figures/total_MAGs_qual_1tig.svg")
