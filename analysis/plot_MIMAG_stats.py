# -*- coding: utf-8 -*-

"""
Created on Mon Sep 02 10:41:58 2024
@author: Anna Cusco
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

### Import data
MIMAG_report = pd.read_csv('/data/Projects/ShanghaiDogs/data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv', \
                    delimiter=',', header=0, index_col=0)
GTDB_qual = pd.read_csv('/data/Projects/ShanghaiDogs/external-data/data/NCBI_genomes_ref/NCBI_genomes_qual_MIMAG_report.csv',sep=',')

MIMAG_report['Quality_det']=MIMAG_report['Quality']
for n in MIMAG_report.index:
    if MIMAG_report.loc[n, 'Nr contigs'] == 1 and MIMAG_report.loc[n, 'Quality'] == 'high-quality':
        MIMAG_report.loc[n, 'Quality_det'] = 'single-contig HQ'

GTDB_qual['Quality_det']=GTDB_qual['Quality']
for n in GTDB_qual.index:
    if GTDB_qual.loc[n, 'Number'] == 1 and GTDB_qual.loc[n, 'Quality'] == 'high-quality':
        GTDB_qual.loc[n, 'Quality_det'] = 'single-contig HQ'

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

### TAXONOMY OF species-level MAGs

species_catalog = MIMAG_report[MIMAG_report['Representative']=='Yes']
species_catalog['ref_new']=species_catalog['GTDBtk fastani Ref']
species_catalog['ref_new'] = species_catalog['ref_new'].fillna('Novel species')

for n in species_catalog.index:
    ref_value = species_catalog.loc[n, 'ref_new']
    if 'GCF' in ref_value:
        species_catalog.loc[n, 'ref_new'] = 'RefSeq reference (GCF)'
    elif 'GCA' in ref_value:
        species_catalog.loc[n, 'ref_new'] = 'GenBank reference (GCA)'

species_catalog['Phylum'] = species_catalog['Classification'].str.extract(r'p__([^;]+)')

### Phylum-level MAG counts by reference genome / novelty
phylum_ref = species_catalog.pivot_table(
    index='Phylum',
    columns='ref_new',
    aggfunc='size',
    fill_value=0
)

phylum_ref['Total']=phylum_ref.sum(axis=1)
phylum_ref=phylum_ref.sort_values(by='Total')
phylum_ref = phylum_ref[['RefSeq reference (GCF)','GenBank reference (GCA)','Novel species']]

# Plot horizontal stacked barplot
# Fig_size
width_mm = 130
height_mm = 70
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

ax = phylum_ref.plot(kind='barh', stacked=True,  figsize=figsize_inch,
                     color=['#a6761d','#e6ab02','#66a61e'],width=0.75)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(title='', loc='lower right',fontsize=10)
ax.set_ylabel('Phylum',fontsize=10)
ax.set_xlabel('N# of species-level MAGs',fontsize=10)

plt.tight_layout()
#plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/analysis/figures/sp_MAG_novelty.svg")

### Phylum-level MAG counts by quality
species_catalog['Quality_det']=species_catalog['Quality']+"_"+species_catalog['MIMAG']
species_catalog['Quality_det']=species_catalog['Quality_det'].str.replace('medium-quality_No','Medium-quality MAGs')
species_catalog['Quality_det']=species_catalog['Quality_det'].str.replace('high-quality_No','High-quality MAGs')
species_catalog['Quality_det']=species_catalog['Quality_det'].str.replace('high-quality_Yes','MIMAG High-quality MAGs')

phylum_qual = species_catalog.pivot_table(
    index='Phylum',
    columns='Quality_det',
    aggfunc='size',
    fill_value=0
)

phylum_qual['Total']=phylum_qual.sum(axis=1)
phylum_qual=phylum_qual.sort_values(by='Total')
phylum_qual = phylum_qual[['MIMAG High-quality MAGs','High-quality MAGs','Medium-quality MAGs']]

# Plot horizontal stacked barplot
# Fig_size
width_mm = 130
height_mm = 70
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

ax = phylum_qual.plot(kind='barh', stacked=True,  figsize=figsize_inch,
                     color=['#1E3F20','#1B9E77','#EDB458'],width=0.75)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(title='', loc='lower right',fontsize=10)
ax.set_ylabel('Phylum',fontsize=10)
ax.set_xlabel('N# of species-level MAGs',fontsize=10)

plt.tight_layout()
#plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/analysis/figures/sp_MAG_quality.svg")
