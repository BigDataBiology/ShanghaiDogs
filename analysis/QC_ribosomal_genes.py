# -*- coding: utf-8 -*-

"""
Created on Tue May 21 11:36:22 2024
@author: Anna Cusco
"""

import os
import glob
import pandas as pd

os.chdir('/data/Projects/ShanghaiDogs/')

# Import files
min_length_rRNA = pd.read_csv('intermediate-outputs/07_ribosomal_genes/ribosomal_gene_ids.tsv',sep='\t')
min_length_rRNA['ribosomal_file']=min_length_rRNA['ribosomal_file'].str.replace('_ribosomal.fa','')

extremes = pd.read_csv('intermediate-outputs/07_ribosomal_genes/ribosomals_extreme_location.csv',index_col=0)
QC_ribosomals = pd.merge(min_length_rRNA,extremes,left_on='ribosomal_file',right_on='Bin ID',how='left').fillna(0)
QC_ribosomals.drop(['ribosomal_file','Bin ID_y','min_id_5s'],inplace=True,axis=1)
QC_ribosomals.columns = ['Bin ID','min_id_16S','mean_id_16S','min_id_23S','16S_extreme','23S_extreme','jug_idx']

MIMAG = pd.read_csv('data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv').fillna(0)
MIMAG['Bin ID']=MIMAG['Bin ID'].str.replace('.fna.gz','')

# Merge reports
MIMAG_QC_rRNA = pd.merge(MIMAG,QC_ribosomals,left_on='Bin ID',right_on='Bin ID')
