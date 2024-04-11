# -*- coding: utf-8 -*-

"""
Created on Tue Feb 20 13:43:05 2024
@author: Anna Cusco
"""

import sys
import pandas as pd
import os

# First merge all the samples qual reports:
# cd /data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool
# head -n 1 D000/D000_qual_report.csv > ALL_bins_qual_report.csv
# awk 'FNR > 1' D0*/D0*_qual_report.csv >> ALL_bins_qual_report.csv

# Import quality report for all bins
qual_report = "/data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool/ALL_bins_qual_report.csv"
qual_report = pd.read_csv(qual_report, delimiter=',',header=0,index_col=0)
qual_report = qual_report[~qual_report['Quality'].str.contains('low_quality')]

# Import list of species-level bins
sp_level = pd.read_csv("/data/Projects/ShanghaiDogs/data/ShanghaiDogsTables/ShanghaiDogsMAGs_ANI95_sp.txt",\
                     header=None,names=["Complete ID"])
sp_level['Representative']='Yes'
qual_report = pd.merge(qual_report,sp_level,how='outer',right_on='Complete ID',left_on='Complete ID')
qual_report.loc[qual_report['Representative'].isna(), 'Ref'] = 'No'

# Import renaming file
renamed_bins = pd.read_csv("/data/Projects/ShanghaiDogs/data/ShanghaiDogsMAGs/mag_meta.tsv.gz",sep='\t')
renamed_bins['Complete ID'] = renamed_bins['Name'] + '_' + renamed_bins['Sample']
renamed_bins = renamed_bins[['Complete ID','Filename']]
qual_report_bins_renamed = pd.merge(qual_report,renamed_bins,left_on='Complete ID',right_on='Complete ID')

# Import and update quality values for Allobaculum stercoricanis (CheckM2 general model, instead of specific - this last is default)
checkm_allobaculum = pd.read_csv("/data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool/Allobaculum_CheckM2_general_model/Allo_quality_report.tsv",\
                                 delimiter='\t',header=0)
qual_report_intermediate = pd.merge(qual_report_bins_renamed,checkm_allobaculum[['Completeness_General','Contamination','Name']],\
                                    how='outer',right_on='Name',left_on='Complete ID')

qual_report_intermediate.loc[qual_report_intermediate['Completeness_General'].isna(), 'Completeness_General'] = qual_report_intermediate['Completeness']
qual_report_intermediate.loc[qual_report_intermediate['Contamination_y'].isna(), 'Contamination_y'] = qual_report_intermediate['Contamination_x']
qual_report_intermediate.drop(['Completeness','Contamination_x','Name'],axis=1,inplace=True)
qual_report_intermediate.rename(columns={'Completeness_General': 'Completeness'}, inplace=True)
qual_report_intermediate.rename(columns={'Contamination_y': 'Contamination'}, inplace=True)

qual_report_intermediate.loc[(qual_report_intermediate['Completeness'] >= 50) & (qual_report_intermediate['Contamination'] < 10), 'Quality'] = 'medium-quality'
qual_report_intermediate.loc[(qual_report_intermediate['Completeness'] >= 90) & (qual_report_intermediate['Contamination'] < 5), 'Quality'] = 'high-quality'

# Import tRNA count for all bins
tRNA_path = "/data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool/tRNA_scan_output"
tRNA_count = pd.DataFrame(columns=['complete_id', 'unique_tRNAs', 'total_tRNAs'])

for folder_name in os.listdir(tRNA_path):
    if os.path.isdir(os.path.join(tRNA_path, folder_name)):
        files_in_folder = os.listdir(os.path.join(tRNA_path, folder_name))
        #print(files_in_folder)
        for file_name in files_in_folder:
            print(file_name)
            bin_path = os.path.join(tRNA_path, folder_name) + "/" + file_name
            tRNA_out = pd.read_csv(bin_path, skiprows=3, delimiter='\t',header=None)
            tRNA_out.columns = ("contig_id","tRNA","Begin","End","Type","Codon","Begin","End","Score","Note")
            tRNA_out['complete_id']=file_name.split('_trna.out')[0]
            aa_count = tRNA_out['Type'].value_counts() # Count each AA Type
            aa_count = aa_count.reset_index()
            unique_aa = len(aa_count['Type'])
            total_aa = aa_count['count'].sum()
            tRNA_count.loc[tRNA_out['complete_id'][0]] = [tRNA_out['complete_id'][0], unique_aa, total_aa]

tRNA_count.drop(['complete_id'],axis=1,inplace=True)

# Create quality report final

qual_report_final=pd.merge(qual_report_intermediate,tRNA_count,left_on='Complete ID',right_index=True)
qual_report_final.columns = ['Classification','GTDBtk fastani Ref','Nr contigs','Quality','Sample','Original ID',
                             '16S','23S','5S','16S partial','23S partial','5S partial','Representative','Bin ID',\
                             'Completeness','Contamination','Unique tRNAs','Total tRNAs']

# Loop through each row in the DataFrame
# Assess if the MAGs follow MIMAG criteria
qual_report_final['MIMAG']='No'

for index, row in qual_report_final.iterrows():
    if row['Quality'] == 'high-quality' and row['16S'] > 0 and row['Unique tRNAs'] > 19:
            print(row['Quality'])
            qual_report_final.loc[index, 'MIMAG'] = 'Yes'

qual_report_final = qual_report_final[['Bin ID','Completeness','Contamination','Quality','Nr contigs','Classification',\
                                       'GTDBtk fastani Ref','16S','23S','5S','16S partial','23S partial','5S partial',\
                                       'Unique tRNAs','Total tRNAs','MIMAG','Sample','Original ID']]

qual_report_final.to_csv("/data/Projects/ShanghaiDogs/data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv",index=False)
