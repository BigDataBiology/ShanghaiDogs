# -*- coding: utf-8 -*-
# date "+Created on %a %d %b %Y %H:%M:%S" >> 1-Beta_diversity.py

"""
Created on Thu 06 Jun 2024 10:44:06
@author: Anna Cusco
"""

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

os.chdir('/data/Projects/ShanghaiDogs/')

## Import input files
# Import abundance table
abd_tab = pd.read_csv('intermediate-outputs/repbin_coverage_rmean.tsv', sep='\t',index_col=0)
text_to_remove = '_SR_to_95_ANI'
abd_tab.columns = [col.replace(text_to_remove, '') for col in abd_tab.columns]

# Import SHD metadata table
metadata = pd.read_csv('data/ShanghaiDogsMetadata/SH_Dog_metadata_red.csv', sep=',', index_col=0)
metadata['Animal_age_simplified'] = metadata['Age_classification']
metadata['Animal_age_simplified'] = metadata['Animal_age_simplified'].str.replace('Early_senior','Senior')
metadata['Animal_age_simplified'] = metadata['Animal_age_simplified'].str.replace('Geriatric','Senior')
metadata['Animal_age_simplified'] = metadata['Animal_age_simplified'].str.replace('Late_senior','Senior')
metadata['Animal_age_simplified'] = metadata['Animal_age_simplified'].str.replace('Mature_adult','Adult')
metadata['Animal_age_simplified'] = metadata['Animal_age_simplified'].str.replace('Juvenile','Young')
metadata['Animal_age_simplified'] = metadata['Animal_age_simplified'].str.replace('Puppy','Young')
metadata['Animal_age_simplified'] = metadata['Animal_age_simplified'].str.replace('Young_adult','Young')

# Import MIMAG table
MIMAG_tab = pd.read_csv('data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv', sep=',')
MIMAG_tab['Bin ID'] = MIMAG_tab['Bin ID'].str.split('.').str[0]

## Add taxonomy information to abundance table

MIMAG_tab[['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']] = \
    MIMAG_tab['Classification'].str.split(';', expand=True)

# Remove the prefixes
MIMAG_tab['domain'] = MIMAG_tab['domain'].str.replace('d__', '')
MIMAG_tab['phylum'] = MIMAG_tab['phylum'].str.replace('p__', '')
MIMAG_tab['class'] = MIMAG_tab['class'].str.replace('c__', '')
MIMAG_tab['order'] = MIMAG_tab['order'].str.replace('o__', '')
MIMAG_tab['family'] = MIMAG_tab['family'].str.replace('f__', '')
MIMAG_tab['genus'] = MIMAG_tab['genus'].str.replace('g__', '')
MIMAG_tab['species'] = MIMAG_tab['species'].str.replace('s__', '')

for n in range(len(MIMAG_tab)):
    if MIMAG_tab['genus'].loc[n] == '':
        MIMAG_tab['genus'].loc[n] = 'unknown ' + MIMAG_tab['family'].loc[n]
        print(MIMAG_tab['genus'].loc[n])

for n in range(len(MIMAG_tab)):
    if MIMAG_tab['species'].loc[n] == '':
        MIMAG_tab['species'].loc[n] = 'unknown ' + MIMAG_tab['genus'].loc[n]
        print(MIMAG_tab['species'].loc[n])

MIMAG_tab = MIMAG_tab.set_index('Bin ID')
abd_tab_tax = pd.merge(abd_tab,MIMAG_tab['species'],left_index=True,right_index=True)
abd_tab_sp = abd_tab_tax.set_index('species')
abd_tab_sp.to_csv('intermediate-outputs/MAGs_tax_profiling/SHDs-MAG_RA_tab.csv')

## Remove low abundant MAGs and plot species-level clustermap

tax_tab = abd_tab_tax.set_index('species')

def rm_low_abd_tax(tax_tab,max):
    """
    Taxonomy/feature table needs samples as columns.
    Remove those tax with a lower max abundance than
    (max) in at least one sample.
    """
    tax_tab['Max'] = tax_tab.max(axis=1)
    tax_tab_filt = tax_tab.loc[tax_tab['Max'].between(max, 1)]
    tax_tab_filt = tax_tab_filt.drop(['Max'], axis=1)
    tax_tab_filt['Sum'] = tax_tab_filt.sum(axis=1)
    tax_tab_filt_ord = tax_tab_filt.sort_values(by='Sum', ascending=False)
    tax_tab_filt_ord = tax_tab_filt_ord.drop(['Sum'], axis=1)
    tax_tab_filt_ord = tax_tab_filt_ord.T
    tax_tab_filt_ord['Other_tax'] = 1 - tax_tab_filt_ord.sum(axis=1)
    tax_tab_tax_ls = tax_tab_filt_ord.T.index.to_list()
    tax_tab_filt_ord = tax_tab_filt_ord.T
    return tax_tab_filt_ord, tax_tab_tax_ls

tax_tab_filt, tax_tab_filt_ls = rm_low_abd_tax(tax_tab,0.05)

# Plot clustermap of the most abundant / prevalent species in SHD

tax_tab_filt = tax_tab_filt.drop(['Other_tax'], axis=0)

mask = (tax_tab_filt == 0)
CM = sns.clustermap(tax_tab_filt, figsize=(12, 8), row_cluster=True, mask = mask, \
                                      cmap='Blues', linewidth=0.3, vmax=0.2, \
                                      yticklabels=True, xticklabels=True)
CM.ax_cbar.set_position((0.02, 0.85, 0.03, 0.13))
plt.show()

## Remove low abundant MAGs and plot genus-level clustermap
abd_tab_tax = pd.merge(abd_tab,MIMAG_tab['genus'],left_index=True,right_index=True)
tax_tab_genus = abd_tab_tax.groupby('genus').sum().reset_index()
tax_tab_genus = tax_tab_genus.set_index('genus')

genus_tab_filt, genus_tab_filt_ls = rm_low_abd_tax(tax_tab_genus,0.02)

## Plot clustermap of the most abundant / prevalent genera in SHD
mask = (genus_tab_filt == 0)
CM = sns.clustermap(genus_tab_filt, figsize=(12, 14), row_cluster=True, mask = mask, \
                                      cmap='Blues', linewidth=0.3, vmax=0.2, \
                                      yticklabels=True, xticklabels=True)
CM.ax_cbar.set_position((0.02, 0.85, 0.03, 0.13))
plt.show()

# PLOT Clustermap with a column colored by specific variable (e.g. Household)
variable = 'Animal_age_simplified'

genus_tab_filt = genus_tab_filt.T
genus_tab_metadata = pd.merge(genus_tab_filt,metadata[variable],left_index=True,right_index=True)
genus_tab_metadata = genus_tab_metadata.reset_index()

# Identify metadata entries exclusive to a single dog & replace unique entries with 'n=1'
mask_unique_var = genus_tab_metadata[variable].value_counts() == 1
unique_var = mask_unique_var[mask_unique_var].index
genus_tab_metadata[variable] = genus_tab_metadata[variable].apply(lambda x: 'n=1' if x in unique_var else x)
genus_tab_metadata.set_index('index', inplace=True)

# Create a color palette with the same number of colors as unique categories in the metadata column
# Create a dictionary where the key is the category and the values are the colors from the palette
sns.set(font_scale=1.2)
categories = genus_tab_metadata[variable].unique()
network_pal = sns.color_palette("tab10", len(categories))
network_dict = dict(zip(categories, network_pal))
networks = genus_tab_metadata[variable]
network_colors = pd.Series(networks).map(network_dict)

#store column names to list
cols = genus_tab_metadata.columns.values.tolist()
cols.remove(variable)

#plot the clustermap
clustermap = sns.clustermap(genus_tab_metadata[cols], row_colors=(network_colors), cmap='Blues', \
                    linewidth=0.3, figsize=(16, 14), vmax=0.5,\
                    yticklabels=True, xticklabels=True)
clustermap.fig.subplots_adjust(right=0.78,left=0.05)
clustermap.ax_cbar.set_position((0.02, .83, .03, .15)) #position color bar
clustermap.ax_heatmap.set_yticklabels([])
clustermap.ax_heatmap.set_yticks([])

# Create legend outside the clustermap
legend_patches = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=label)
                  for label, color in network_dict.items()]
plt.legend(handles=legend_patches, title=variable, bbox_to_anchor=(32, 1), loc='upper right')

fig_path='intermediate-outputs/MAGs_tax_profiling/figs/clustermap_'+variable+'.svg'
clustermap.fig.savefig(fig_path, format='svg')
plt.show()