# -*- coding: utf-8 -*-

"""
Created on Mon Sep 02 10:41:58 2024
@author: Anna Cusco
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import squarify

plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

### Import data
MIMAG_report = pd.read_csv('/data/Projects/ShanghaiDogs/data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv', \
                    delimiter=',', header=0, index_col=0)
GTDB_qual = pd.read_csv('/data/Projects/ShanghaiDogs/external-data/data/NCBI_genomes_ref/NCBI_genomes_qual_MIMAG_report.csv',sep=',')

MIMAG_report['Quality_det']=MIMAG_report['Quality']
for n in MIMAG_report.index:
    if MIMAG_report.loc[n, 'Nr contigs'] == 1 and MIMAG_report.loc[n, 'Quality'] == 'high-quality':
        MIMAG_report.loc[n, 'Quality_det'] = 'single-contig HQ'

MIMAG_report['Quality_MIMAG']=MIMAG_report['Quality']
for n in MIMAG_report.index:
    if MIMAG_report.loc[n, 'MIMAG'] == 'Yes' and MIMAG_report.loc[n, 'Quality'] == 'high-quality':
        MIMAG_report.loc[n, 'Quality_MIMAG'] = 'MIMAG high-quality'

GTDB_qual['Quality_det']=GTDB_qual['Quality']
for n in GTDB_qual.index:
    if GTDB_qual.loc[n, 'Number'] == 1 and GTDB_qual.loc[n, 'Quality'] == 'high-quality':
        GTDB_qual.loc[n, 'Quality_det'] = 'single-contig HQ'

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
species_catalog['Family'] = species_catalog['Classification'].str.extract(r'f__([^;]+)')
species_catalog['Genus'] = species_catalog['Classification'].str.extract(r'g__([^;]+)')
species_catalog['Species'] = species_catalog['Classification'].str.extract(r's__([^;]+)')

### HQ_mq MAGs table
quality_table = MIMAG_report.pivot_table(
    index='Sample',
    columns='Quality',
    aggfunc='size',
    fill_value=0)

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
plt.savefig("/data/Projects/ShanghaiDogs/intermediate-outputs/figures/total_MAGs_qual.svg")

### HQ_mq 1-contig
quality_table_det = MIMAG_report.pivot_table(
    index='Sample',
    columns='Quality_det',
    aggfunc='size',
    fill_value=0)

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
#plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/intermediate-outputs/figures/total_MAGs_qual_1tig.svg")

### MIMAG high-quality
quality_table_MIMAG = MIMAG_report.pivot_table(
    index='Sample',
    columns='Quality_MIMAG',
    aggfunc='size',
    fill_value=0)

quality_table_MIMAG['Total']=quality_table_MIMAG['high-quality']+quality_table_MIMAG['medium-quality']+quality_table_MIMAG['MIMAG high-quality']
quality_table_MIMAG = quality_table_MIMAG.sort_values('Total',ascending=False)

# Plot a stacked barplot
width_mm = 110
height_mm = 60
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

# Barplot
ax = (quality_table_MIMAG[['MIMAG high-quality','high-quality', 'medium-quality']].plot
      (kind='bar', stacked=True,color=['#1B9E77','#a6761d','#EDB458'],width=0.75,
       figsize=figsize_inch))

ax.set_ylabel('N# of MAGs',fontsize=10)
ax.set_xlabel('Samples',fontsize=10)
ax.set_yticklabels(map(int, ax.get_yticks()), fontsize=10)
ax.set_xticklabels([])
ax.tick_params(axis='x', bottom=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Custom legend text
legend_labels = ['MIMAG High-quality MAGs', 'High-quality MAGs', 'Medium-quality MAGs']
ax.legend(labels=legend_labels, title='', loc='upper right',fontsize=10)

# Save figure
plt.tight_layout()
#plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/intermediate-outputs/figures/total_MAGs_qual_MIMAG.svg")

### Most prevalent species in dog cohort
# Create 'species' column
MIMAG_report['species'] = MIMAG_report['Classification'].str.extract(r's__([^;]+)')
print(MIMAG_report[['Classification', 'species']].head())
MIMAG_report_species = MIMAG_report.dropna(subset=['species'])

# Create a pivot table with counts based on 'Family' and 'Quality'
prevalent_sp = MIMAG_report_species.pivot_table(
    index='species',
    columns='Quality_det',
    values='Classification',
    aggfunc='count').fillna(0)

prevalent_sp['Total']=prevalent_sp['high-quality']+prevalent_sp['medium-quality']+prevalent_sp['single-contig HQ']
prevalent_sp = prevalent_sp.sort_values(by='Total', ascending=False).head(25)
prevalent_sp_perc = prevalent_sp.div(prevalent_sp['Total'], axis=0) * 100
prevalent_sp_perc = prevalent_sp_perc.sort_values(by='medium-quality',ascending=False)

prevalent_sp_perc_phylum = pd.merge(prevalent_sp_perc,species_catalog[['Species','Phylum']],left_index=True,right_on='Species')
order = ['Bacillota_A',
 'Bacteroidota',
 'Bacillota',
 'Fusobacteriota',
 'Pseudomonadota',
 'Actinomycetota',
 'Bacillota_C',
 'Campylobacterota',
 'Bacillota_B',
 'Desulfobacterota',
 'Deferribacterota']

prevalent_sp_perc_phylum['Phylum'] = pd.Categorical(
    prevalent_sp_perc_phylum['Phylum'],
    categories=order,
    ordered=True)

prevalent_sp_perc_phylum = prevalent_sp_perc_phylum.sort_values(
    by=['Phylum', 'Species'],
    ascending=[False, False])

prevalent_sp_perc_phylum = prevalent_sp_perc_phylum.set_index('Species')

# Plotting
width_mm = 95
height_mm = 120
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)
prevalent_sp_perc_phylum[['single-contig HQ','high-quality', 'medium-quality']].plot(
    kind='barh',
    stacked=True,
    color=['#1E3F20','#1B9E77', '#EDB458'],
    ax=ax,
    width=0.7)

ax.set_xlabel('',fontsize=10)
ax.set_xticks([0, 50, 100])
ax.set_xticklabels(['0%','50%','100%'])
ax.set_ylabel('')
ax.set_xlim(0,100)

ax.get_legend().remove()
plt.tight_layout()
#plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/intermediate-outputs/figures/prevalent_sp-MAGs_qual_30.svg")

### Phylum-level MAG counts by reference genome / novelty
phylum_ref = species_catalog.pivot_table(
    index='Phylum',
    columns='ref_new',
    aggfunc='size',
    fill_value=0)

phylum_ref['Total']=phylum_ref.sum(axis=1)
phylum_ref=phylum_ref.sort_values(by='Total',ascending=True)
phylum_ref = phylum_ref[['RefSeq reference (GCF)','GenBank reference (GCA)','Novel species']]
phylum_ref = phylum_ref.tail(8)
phylum_ref.columns = ['RefSeq reference', 'GenBank reference', 'Novel species']

# Plot horizontal stacked barplot
width_mm = 110
height_mm = 65
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

ax = phylum_ref.plot(kind='barh', stacked=True,  figsize=figsize_inch,
                     color=['#a6761d','#e6ab02','#66a61e'],width=0.75)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(title='', loc='lower right',fontsize=10)
ax.set_ylabel('',fontsize=10)
ax.set_xlabel('N# of species-level MAGs',fontsize=10)

plt.tight_layout()
#plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/intermediate-outputs/figures/sp_MAG_novelty_red.svg")

### NOVEL SPECIES PLOTS
species_catalog_novel = species_catalog.query('ref_new == "Novel species"')
species_catalog_novel = species_catalog_novel[['Family','Genus']]

novel_df = species_catalog_novel.pivot_table(
    index='Family',
    aggfunc='count',
    fill_value=0)

novel_df = novel_df.sort_values(by='Genus',ascending=False)

# Create other_tax category
other_tax = novel_df[novel_df['Genus'] < 5]['Genus'].sum()
other_row = pd.DataFrame({'Family': ['Other tax (<5 species)'], 'Genus': [other_tax]})

# Filter out rows with counts less than 4
novel_df = novel_df.reset_index()
novel_df_filtered = novel_df[novel_df['Genus'] >= 5]

# Append the "Other tax" row
final_df = pd.concat([novel_df_filtered, other_row], ignore_index=True)

## DONUT PLOT
width_mm = 110
height_mm = 50
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)
ax.clear()

genus_counts=list(final_df['Genus'])
family_labels=list(final_df['Family'])

# Create a pieplot
plt.pie(genus_counts,labels=family_labels,
        colors=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#D3D3D3'],
        wedgeprops = { 'linewidth' : 3, 'edgecolor' : 'white' },
        textprops={'fontsize': 10})

# add a circle at the center to transform it in a donut chart
my_circle=plt.Circle( (0,0), 0.5, color='white')
p=plt.gcf()
p.gca().add_artist(my_circle)

plt.tight_layout()
#plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/intermediate-outputs/figures/novel_tax_donutplot.svg")

## TREE PLOT
# Set up the figure
width_mm = 160
height_mm = 50
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)
ax.clear()

# Draw the treemap
squarify.plot(sizes=genus_counts,
              label=family_labels,
              ax=ax,
              color=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#D3D3D3'],
              text_kwargs={'fontsize':9})

# Customize the plot
ax.set_axis_off()  # Hide the axes

# Display the plot
plt.tight_layout()
#plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/intermediate-outputs/figures/novel_tax_squarify.svg")

## GENUS PIECHARTS FOR NOVEL SPECIES
all_df_genus = species_catalog.groupby(['Genus','ref_new'])['Representative'].count().reset_index()
all_df_genus_wide = (all_df_genus.pivot_table(index='Genus',
                                             columns='ref_new',
                                             values='Representative')
                     .fillna(0))

all_df_genus_wide['Total'] = all_df_genus_wide.sum(axis=1)
all_df_genus_wide_filt = all_df_genus_wide[all_df_genus_wide['Novel species'] > 3]
all_df_genus_wide_filt.columns = ['GenBank','Novel','RefSeq','Total']

## PIE CHART
genus = 'Dysosmobacter'  # CAG-269, Collinsella, Blautia_A, Dysosmobacter
custom_palette = {'GenBank': '#a6761d', 'Novel': '#66a61e', 'RefSeq': '#e6ab02'}

# Size and dimensions of the pie chart
size_factor = all_df_genus_wide_filt.loc[genus, 'Total']
width_mm = 3 * size_factor
height_mm = 3 * size_factor
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

# Transpose the data for the genus, excluding 'Total'
all_df_genus_wide_filt_T = all_df_genus_wide_filt.drop(['Total'], axis=1).T

# Prepare the data for plotting
genus_counts = list(all_df_genus_wide_filt_T[genus])

# Create the pie chart
fig, ax = plt.subplots(figsize=figsize_inch)
ax.clear()
plt.pie(
    genus_counts,
    labels=None,
    colors=[custom_palette[col] for col in all_df_genus_wide_filt_T.index],
    wedgeprops={'linewidth': 2, 'edgecolor': 'white'})

#plt.show()
out_path = "/data/Projects/ShanghaiDogs/intermediate-outputs/figures/novel_sp_"+genus+".svg"
plt.savefig(out_path)
