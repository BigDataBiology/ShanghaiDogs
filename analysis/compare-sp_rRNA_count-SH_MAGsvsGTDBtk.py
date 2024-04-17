# -*- coding: utf-8 -*-

"""
Created on Fri Apr 12 12:07:54 2024
@author: Anna Cusco
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

os.chdir('/data/Projects/ShanghaiDogs/')

# Import quality reports
SHD_qual = pd.read_csv('data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv',sep=',')
GTDB_qual = pd.read_csv('external-data/data/NCBI_genomes_ref/NCBI_genomes_qual_MIMAG_report.csv',sep=',')

# Filter 'MIMAG' genomes from Shanghai dogs
SHD_qual_HQ = SHD_qual[SHD_qual['Quality']=='high-quality']
count_sp = SHD_qual_HQ.groupby('Classification').size().reset_index()

# Create Genus and Species column
count_sp['Genus'] = count_sp['Classification'].str.replace(r'^.*g__', '', regex=True)
count_sp['Genus'] = count_sp['Genus'].str.replace(r';s__.*', '', regex=True)
count_sp['Species'] = count_sp['Classification'].str.replace(r'^.*s__', '', regex=True)
count_sp['Species'] = count_sp['Species'].replace('', 'unclassified')

abd_sp = count_sp[(count_sp[0] >= 20) & (count_sp['Species'] != 'unclassified')]
abd_sp_ls = list(abd_sp['Species'])

# Stripplot of certain species

for sp in abd_sp_ls:
    filt = SHD_qual_HQ[SHD_qual_HQ['Classification'].str.contains(sp)]
    if not filt.empty:
        ref_id = filt['GTDBtk fastani Ref'].iloc[0]
        ref_values = GTDB_qual[GTDB_qual['Name'].str.contains(ref_id)]
        fig, ax = plt.subplots()
        ax.clear()
        sns.stripplot(data=[filt['16S total'],filt['23S total'],filt['5S total']], palette="Dark2",ax=ax,alpha=0.5)
        extra_values = [ref_values['16S total'].iloc[0], ref_values['23S total'].iloc[0], ref_values['5S total'].iloc[0]]
        sns.scatterplot(x=[0, 1, 2], y=extra_values, color='black', ax=ax, zorder=3, marker='X')
        plt.ylabel('Nr ribosomal genes')
        plt.xticks(ticks=[0, 1, 2],labels=['16S rRNA', '23S rRNA','5S rRNA'])
        plt.title(sp)
        sns.despine(fig, trim=False)
        #plt.show()
        out_path='analysis/figures/rRNA_count_species/rRNAs_'+sp+'.png'
        plt.savefig(out_path)
    else:
        print(sp + ' is not in the GTDB_qual table')
