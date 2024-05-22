# -*- coding: utf-8 -*-

"""
Created on Mon May 17 10:12:01 2024
@author: Anna Cusco
"""

import os
import glob
import pandas as pd

os.chdir('/data/Projects/ShanghaiDogs/')

# Import contig size (Flye - it might be slightly different after polypolish but should work)
assembly_path = '/data/Projects/ShanghaiDogs/data/ShanghaiDogsAssemblies/'
assembly_info = []

for file in glob.glob(os.path.join(assembly_path, '*.txt')):
    try:
        info = pd.read_csv(file, delimiter='\t')
        info = info[['#seq_name','length']]
        info['sample']=os.path.basename(file).split('.')[0]
        info['sample']=info['sample'].str.split('_').str[0]
        assembly_info.append(info)
        #print(info.head())
    except Exception as e:
        print(f'Error reading {file}: {e}')

assembly_all = pd.concat(assembly_info, ignore_index=True)

### Import subset_assemblies info
assembly_path_sub = '/data/anna/animal_metagenome/long-mg-dog/01_assembly/subset_assemblies_20G/*/assembly_info.txt'
assembly_info = []

for file in glob.glob(assembly_path_sub):
    try:
        sample_name = file.split('/')[7].split('_')[0] + '_10G_20G'
        info = pd.read_csv(file, delimiter='\t', usecols=['#seq_name', 'length'])
        info['sample'] = sample_name
        info['complete_id'] = info['sample'] + '_' + info['#seq_name']
        assembly_info.append(info)
        print(info.head())
    except Exception as e:
        print(f'Error reading {file}: {e}')

assembly_sub = pd.concat(assembly_info, ignore_index=True)
assembly_sub_max_length = assembly_sub.loc[assembly_sub.groupby('complete_id')['length'].idxmax()]
assembly_sub_final = assembly_sub_max_length.drop_duplicates(subset='complete_id', keep='last')
assembly_sub_final.drop(['complete_id'],inplace=True,axis=1)

assembly_merged = pd.concat([assembly_all,assembly_sub_final],ignore_index=True)

# Import rRNA ribosomal files
rRNA_path = '/data/Projects/ShanghaiDogs/intermediate-outputs/07_ribosomal_genes/barrnap_out/'
rRNAs_out = []
for folder_name in os.listdir(rRNA_path):
    file_path = rRNA_path + folder_name + '/' + folder_name + '_ALL.txt'
    print(file_path)
    try:
        df = pd.read_csv(file_path, delimiter='\t',header=None)
        df['sample'] = df[0].apply(lambda x: x.split('_')[0] + '_10G_20G' if '10G_20G' in x else x.split('_')[0])
        print(df['sample'])
        df['contig'] = df[1].str.split('_polypolish').str[0]
        df = df[['sample','contig',0,4,5,9]]
        df.columns = ['sample','contig','bin_id','start','end','annotation']
        rRNAs_out.append(df)
        print(df.head())
    except Exception as e:
        print(f'Error reading {file_path}: {e}')

rRNAs_all = pd.concat(rRNAs_out, ignore_index=True)

# Import jug_index
jug_idx = pd.read_csv('intermediate-outputs/07_ribosomal_genes/jug_idxs.csv')
jug_idx['Bin ID'] = jug_idx['ribosomal_files'].str.split('/').str[-1]
jug_idx['Bin ID'] = jug_idx['Bin ID'].str.replace('_ribosomal.fa','')
jug_idx.drop(['ribosomal_files'],inplace=True,axis=1)

# Merge the two datasets
merged = pd.merge(rRNAs_all,assembly_merged,left_on=['sample','contig'],right_on=['sample','#seq_name'],how='left')
merged['Bin ID']=merged['bin_id'].str.replace('_barrnap.txt','')
merged = merged.sort_values(by='bin_id')
merged = merged[~merged['annotation'].str.contains('5S')] #ignore 5S rRNA genes

# Calculate how close are the ribosomal from the 'extremes' of the contig
merged['dist_end']=merged['length']-merged['end']
merged['dist_start']=merged['start']

# Removal of D033_10G_20G, D041_10G_20G and D049_10G_20G
rm_samples=['D033_10G_20G','D041_10G_20G','D049_10G_20G']
merged_filt=merged[~merged['sample'].isin(rm_samples)]

# keep only ribosomal info for partial, and rRNA genes close to the extremes
partial_rRNAs = merged_filt[merged_filt['annotation'].str.contains('partial')]
filtered_16S = merged_filt[(merged_filt['annotation'].str.contains('16S')) & \
                           ((merged_filt['dist_start'] < 1600) | (merged_filt['dist_end']< 1600))]
filtered_23S = merged_filt[(merged_filt['annotation'].str.contains('23S')) & \
                           ((merged_filt['dist_start'] < 3300) | (merged_filt['dist_end'] < 3300))]

merged_filt_extremes = pd.concat([partial_rRNAs,filtered_16S, filtered_23S])
merged_filt_extremes = merged_filt_extremes.drop_duplicates()

# Summarize by Bin ID

counts_16S = merged_filt_extremes[merged_filt_extremes['annotation'].str.contains('16S')].groupby(['Bin ID']).size().reset_index(name='16S_extreme')
counts_23S = merged_filt_extremes[merged_filt_extremes['annotation'].str.contains('23S')].groupby(['Bin ID']).size().reset_index(name='23S_extreme')
extremes_counts = pd.merge(counts_16S, counts_23S, on='Bin ID', how='outer').fillna(0)
extremes_counts = pd.merge(extremes_counts,jug_idx,left_on='Bin ID',right_on='Bin ID')

extremes_counts.to_csv('intermediate-outputs/07_ribosomal_genes/ribosomals_extreme_location.csv')