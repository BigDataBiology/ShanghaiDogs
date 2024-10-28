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

### MIMAG high-quality

quality_table_MIMAG = MIMAG_report.pivot_table(
    index='Sample',
    columns='Quality_MIMAG',
    aggfunc='size',
    fill_value=0
)

quality_table_MIMAG['Total']=quality_table_MIMAG['high-quality']+quality_table_MIMAG['medium-quality']+quality_table_MIMAG['MIMAG high-quality']
quality_table_MIMAG = quality_table_MIMAG.sort_values('Total',ascending=False)

# Plot a stacked barplot
# Fig_size
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
# plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/analysis/figures/total_MAGs_qual_MIMAG.svg")

### Most prevalent species in dog cohort
# Create 'species' column
MIMAG_report['species'] = MIMAG_report['Classification'].str.extract(r's__([^;]+)')
print(MIMAG_report[['Classification', 'species']].head())
MIMAG_report_species = MIMAG_report.dropna(subset=['species'])

# Create a pivot table with counts based on 'Family' and 'Quality'
prevalent_sp = MIMAG_report_species.pivot_table(
    index='species',
    columns='Quality',
    values='Classification',
    aggfunc='count').fillna(0)

prevalent_sp['Total']=prevalent_sp['high-quality']+prevalent_sp['medium-quality']
prevalent_sp = prevalent_sp.sort_values(['Total'],ascending=False)

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

### CONTIGUITY:  MAGs vs REF (GCA & GCF)
# species_catalog_HQ=species_catalog.query('Quality == "high-quality" and ref_new != "Novel species"')
order = ['REF_RefSeq reference (GCF)','MAG_RefSeq reference (GCF)','REF_GenBank reference (GCA)','MAG_GenBank reference (GCA)']
color_palette = ['#a6761d', '#1b9e77', '#e6ab02', '#1b9e77']

### Contiguity by Ref quality vs MAG
SHD_contiguity_df = species_catalog[['Nr contigs','ref_new','GTDBtk fastani Ref']].reset_index()
ALL_contiguity_df = pd.merge(SHD_contiguity_df,GTDB_qual[['Name','Number']],right_on='Name',left_on='GTDBtk fastani Ref')
ALL_contiguity = ALL_contiguity_df [['Bin ID','Nr contigs','Number','ref_new']]
ALL_contiguity.columns = ['Bin ID','MAG','REF','ref_new']

ALL_contiguity_melted = pd.melt(ALL_contiguity, id_vars=['Bin ID','ref_new'],
                    value_vars=['MAG', 'REF'],
                    var_name='Genome', value_name='Count')
ALL_contiguity_melted['category']=ALL_contiguity_melted['Genome']+'_'+ALL_contiguity_melted['ref_new']
ALL_contiguity_melted['category'] = pd.Categorical(ALL_contiguity_melted['category'], categories=order, ordered=True)
ALL_contiguity_melted = ALL_contiguity_melted.sort_values('category')
ALL_contiguity_melted['log count'] = np.log10(ALL_contiguity_melted['Count'])

# Plot boxplot contiguity GCA vs GCF
width_mm = 55
height_mm = 40
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)
ax.clear()
sns.boxplot(data=ALL_contiguity_melted,
              x='category', y='log count',
              palette=['#a6761d', '#1b9e77', '#e6ab02','#1b9e77'],
              ax=ax,
               width=0.8,
               linewidth=1,
               flierprops={
                   'marker': 'd',  # Shape of outliers
                   'color': 'gray',  # Color of outliers
                   'markersize': 2.5,  # Size of outliers
                   'linestyle': 'none'  # No connecting line for outliers
                })
ax.set_ylabel('Log10 counts')
ax.set_ylim(-0.2,3)
ax.set_xlabel('')
ax.set_xticklabels([])
ax.tick_params(axis='x', bottom=False)
sns.despine(fig, trim=False)

plt.tight_layout()
# plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/analysis/figures/sp_MAG_vs_ref_contiguity_boxplot_log_mq-HQ.svg")

### RIBOSOMALS:  MAGs vs REF (GCA & GCF)
order = ['REF_RefSeq reference (GCF)','MAG_RefSeq reference (GCF)','REF_GenBank reference (GCA)','MAG_GenBank reference (GCA)']
color_palette = ['#a6761d', '#1b9e77', '#e6ab02', '#1b9e77']

SHD_ribosomal_df = species_catalog[['16S rRNA','ref_new','GTDBtk fastani Ref']].reset_index()
ALL_ribosomal_df = pd.merge(SHD_ribosomal_df,GTDB_qual[['Name','16S rRNA']],right_on='Name',left_on='GTDBtk fastani Ref')
ALL_ribosomal = ALL_ribosomal_df [['Bin ID','16S rRNA_x','16S rRNA_y','ref_new']]
ALL_ribosomal.columns = ['Bin ID','MAG','REF','ref_new']

ALL_ribosomal_melted = pd.melt(ALL_ribosomal, id_vars=['Bin ID','ref_new'],
                    value_vars=['MAG', 'REF'],
                    var_name='Genome', value_name='Count')
ALL_ribosomal_melted['category']=ALL_ribosomal_melted['Genome']+'_'+ALL_ribosomal_melted['ref_new']
ALL_ribosomal_melted['category'] = pd.Categorical(ALL_ribosomal_melted['category'], categories=order, ordered=True)
ALL_ribosomal_melted = ALL_ribosomal_melted.sort_values('category')

# Plot boxplot 16S rRNA GCA vs GCF
width_mm = 55
height_mm = 40
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)
ax.clear()
sns.boxplot(data=ALL_ribosomal_melted,
            x='category',y='Count',
            palette=['#a6761d', '#1b9e77', '#e6ab02','#1b9e77'],
            ax=ax,
            width=0.7,
            linewidth=1,
            flierprops={
                'marker': 'd',  # Shape of outliers
                'color': 'gray',  # Color of outliers
                'markersize': 2,  # Size of outliers
                'linestyle': 'none'  # No connecting line for outliers
    })
ax.set_ylabel('')
ax.set_xlabel('')
ax.set_xticklabels([])
ax.tick_params(axis='x', bottom=False)
sns.despine(fig, trim=False)

plt.tight_layout()
# plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/analysis/figures/sp_MAG_vs_ref_16S_boxplot.svg")

### Plot histograms for 16S rRNA count
categories = ALL_ribosomal_melted['category'].unique()
num_categories = len(categories)

# Create subplots
width_mm = 120
height_mm = 60
figsize_inch_hist = (width_mm / 25.4, height_mm / 25.4)

plt.clf()
fig, axes = plt.subplots(1, 2, figsize=figsize_inch_hist, sharey=True, sharex=True)
bin_edges = np.linspace(0, 18, 8) # Define common bin edges

# Plot the first two categories on the first subplot
subset1 = ALL_ribosomal_melted[ALL_ribosomal_melted['category'] == categories[0]]
subset2 = ALL_ribosomal_melted[ALL_ribosomal_melted['category'] == categories[1]]
axes[0].hist(subset1['Count'], bins=bin_edges, alpha=0.7, color='#a6761d', edgecolor='white',linewidth=1.2)
axes[0].hist(subset2['Count'], bins=bin_edges, alpha=0.5, color='#1b9e77', edgecolor='white',linewidth=1.2)
axes[0].set_title('')
axes[0].set_xlabel('16S rRNA counts',fontsize=10)
axes[0].set_ylabel('Frequency',fontsize=10)
axes[0].set_xlim(-0.2,17.7)

# Plot the next two categories on the second subplot
subset3 = ALL_ribosomal_melted[ALL_ribosomal_melted['category'] == categories[2]]
subset4 = ALL_ribosomal_melted[ALL_ribosomal_melted['category'] == categories[3]]
axes[1].hist(subset3['Count'], bins=bin_edges, alpha=0.7, color='#e6ab02', edgecolor='white',linewidth=1.2)
axes[1].hist(subset4['Count'], bins=bin_edges, alpha=0.5, color='#1b9e77', edgecolor='white',linewidth=1.2)
axes[1].set_title('')
axes[1].set_xlabel('16S rRNA counts',fontsize=10)
#axes[1].set_ylabel('Frequency',fontsize=10)
axes[1].set_xlim(-0.2,17.7)

sns.despine(fig, trim=False)
plt.tight_layout()
#plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/analysis/figures/sp_MAG_vs_ref_16S_histograms_v2.svg")

### tRNAs: MAGs vs REF (GCA & GCF)
order = ['REF_RefSeq reference (GCF)','MAG_RefSeq reference (GCF)','REF_GenBank reference (GCA)','MAG_GenBank reference (GCA)']
color_palette = ['#a6761d', '#1b9e77', '#e6ab02', '#1b9e77']

SHD_tRNAs_df = species_catalog[['Unique tRNAs','ref_new','GTDBtk fastani Ref']].reset_index()
ALL_tRNAs_df = pd.merge(SHD_tRNAs_df,GTDB_qual[['Name','Unique tRNAs']],right_on='Name',left_on='GTDBtk fastani Ref')
ALL_tRNAs = ALL_tRNAs_df [['Bin ID','Unique tRNAs_x','Unique tRNAs_y','ref_new']]
ALL_tRNAs.columns = ['Bin ID','MAG','REF','ref_new']

ALL_tRNAs_melted = pd.melt(ALL_tRNAs, id_vars=['Bin ID','ref_new'],
                    value_vars=['MAG', 'REF'],
                    var_name='Genome', value_name='Count')
ALL_tRNAs_melted['category']=ALL_tRNAs_melted['Genome']+'_'+ALL_tRNAs_melted['ref_new']
ALL_tRNAs_melted['category'] = pd.Categorical(ALL_tRNAs_melted['category'], categories=order, ordered=True)
ALL_tRNAs_melted = ALL_tRNAs_melted.sort_values('category')

# Plot boxplot tRNAs GCA vs GCF
width_mm = 65
height_mm = 60
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)
ax.clear()
sns.boxplot(data=ALL_tRNAs_melted,
               x='category',y='Count',
               palette=color_palette,
               ax=ax,
               width=0.7,
               linewidth=1,
               flierprops={
                   'marker': 'd',  # Shape of outliers
                   'color': 'gray',  # Color of outliers
                   'markersize': 2,  # Size of outliers
                   'linestyle': 'none'  # No connecting line for outliers
    })
ax.set_ylabel('')
ax.set_xlabel('')
ax.set_xticklabels([])
ax.set_ylim(0,30)
ax.tick_params(axis='x', bottom=False)
sns.despine(fig, trim=False)

plt.tight_layout()
# plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/analysis/figures/sp_MAG_vs_ref_tRNAs_boxplot.svg")

### SWARMPLOTS
fig, ax = plt.subplots(figsize=figsize_inch)
ax.clear()
df=ALL_ribosomal_melted

sns.swarmplot(data=df,
              x='category', y='Count',
              ax=ax,
              hue='category',
              palette=color_palette,
              size=2.5,             # Size of the dots
              alpha=0.5)            # Transparency of the dots

ax.legend_.remove()
ax.set_ylabel('')
ax.set_xlabel('')
ax.set_xticklabels([])
#ax.set_ylim(0, 30)
ax.tick_params(axis='x', bottom=False)
sns.despine(fig, trim=False)

plt.tight_layout()
#plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/analysis/figures/sp_MAG_vs_ref_ribosomal_swarmplot.svg")

### STRIPPLOT
fig, ax = plt.subplots(figsize=figsize_inch)
ax.clear()

df=ALL_tRNAs_melted
categories = df['category'].unique()

# Adding individual data points using stripplot with increased jitter
sns.stripplot(data=df,
              x='category', y='Count',
              ax=ax,
              hue='category',
              palette=color_palette,
              size=3,
              alpha=0.3,  # Transparency of the dots
              jitter=0.2,  # Increase jitter to spread out the dots
              order=categories)

sns.boxplot(data=df,
            x='category', y='Count',
            ax=ax,
            width=0.8,
            palette=[(1, 1, 1, 0) for _ in color_palette], # Transparent boxplot
            linewidth=1,
            fliersize=0) # Hide the outlier markers

ax.legend_.remove()
ax.set_title('Number of tRNAs')
ax.set_ylabel('Counts')
ax.set_xlabel('')
ax.set_xticklabels([])
ax.set_ylim(0, 30)
ax.tick_params(axis='x', bottom=False)
sns.despine(fig, trim=False)

plt.tight_layout()
#plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/analysis/figures/sp_MAG_vs_ref_tRNAs_boxplot-strip.svg")

### NOVEL SPECIES PLOT
species_catalog_novel = species_catalog.query('ref_new == "Novel species"')
species_catalog_novel = species_catalog_novel[['Family','Genus']]

novel_df = species_catalog_novel.pivot_table(
    index='Family',
    aggfunc='count',
    fill_value=0
)

novel_df = novel_df.sort_values(by='Genus',ascending=False)

# Create other_tax category
other_tax = novel_df[novel_df['Genus'] < 5]['Genus'].sum()
other_row = pd.DataFrame({'Family': ['Other tax (<5 species)'], 'Genus': [other_tax]})

# Filter out rows with counts less than 4
novel_df = novel_df.reset_index()
novel_df_filtered = novel_df[novel_df['Genus'] >= 5]

# Append the "Other tax" row
final_df = pd.concat([novel_df_filtered, other_row], ignore_index=True)

# Donut plot
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
my_circle=plt.Circle( (0,0), 0.7, color='white')
p=plt.gcf()
p.gca().add_artist(my_circle)

plt.tight_layout()
#plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/analysis/figures/novel_tax_donutplot.svg")

### Tree plot
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
plt.savefig("/data/Projects/ShanghaiDogs/analysis/figures/novel_tax_squarify.svg")