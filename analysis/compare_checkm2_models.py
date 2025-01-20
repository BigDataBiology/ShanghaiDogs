import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import numpy as np

os.chdir('/data/Projects/ShanghaiDogs')

# Import input files
mimag_report = pd.read_csv('data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv',
                           sep=',', index_col=0)
mimag_report.index = mimag_report.index.str.replace('.fna.gz','')

checkm_polca = pd.read_csv('intermediate-outputs/polish_evaluation_mags/'
                           'polca_bins_checkm_all/quality_report.tsv', sep='\t', index_col=0)
checkm_polca.index = checkm_polca.index.str.replace('.fna','')

checkm_poly = pd.read_csv('intermediate-outputs/polish_evaluation_mags/'
                           'polypolish_bins_checkm_all/quality_report.tsv', sep='\t', index_col=0)
checkm_medaka = pd.read_csv('intermediate-outputs/polish_evaluation_mags/'
                           'medaka_bins_checkm_all/quality_report.tsv', sep='\t', index_col=0)
checkm_flye = pd.read_csv('intermediate-outputs/polish_evaluation_mags/'
                           'flye_bins_checkm_all/quality_report.tsv', sep='\t', index_col=0)

polish_eval = pd.read_csv("intermediate-outputs/polishing_evaluation/results/polishing_coverM_out_LR.csv", \
                  delimiter=',',header=0, index_col=0)

# Create quality_df to check evolution of completeness and contamination across polishing steps
quality_df = pd.merge(mimag_report[['Completeness']], checkm_polca[['Completeness_General', 'Completeness_Specific',
                                                                    'Completeness_Model_Used', 'Contamination']],
                      left_index=True, right_index=True)
quality_df.columns = ['MIMAG-table_compl','Polca-compl_general','Polca-compl_specific','Polca-model','Polca-contam']

quality_df = pd.merge(quality_df,checkm_poly[['Completeness_General', 'Completeness_Specific',
                                              'Completeness_Model_Used', 'Contamination']],
                      left_index=True, right_index=True)
quality_df.columns = ['MIMAG-table_compl', 'Polca-compl_general', 'Polca-compl_specific', 'Polca-model', 'Polca-contam',
                      'Poly-compl_general', 'Poly-compl_specific', 'Poly-model', 'Poly-contam' ]

quality_df = pd.merge(quality_df,checkm_medaka[['Completeness_General', 'Completeness_Specific',
                                                'Completeness_Model_Used', 'Contamination']],
                      left_index=True, right_index=True)
quality_df.columns = ['MIMAG-table_compl', 'Polca-compl_general', 'Polca-compl_specific', 'Polca-model', 'Polca-contam',
                      'Poly-compl_general', 'Poly-compl_specific', 'Poly-model', 'Poly-contam', 'Medaka-compl_general',
                      'Medaka-compl_specific', 'Medaka-model', 'Medaka-contam']

quality_df = pd.merge(quality_df,checkm_flye[['Completeness_General', 'Completeness_Specific',
                                              'Completeness_Model_Used', 'Contamination']],
                      left_index=True, right_index=True)
quality_df.columns = ['MIMAG-table_compl', 'Polca-compl_general', 'Polca-compl_specific', 'Polca-model', 'Polca-contam',
                      'Poly-compl_general', 'Poly-compl_specific', 'Poly-model', 'Poly-contam', 'Medaka-compl_general',
                      'Medaka-compl_specific', 'Medaka-model', 'Medaka-contam', 'Flye-compl_general','Flye-compl_specific',
                      'Flye-model', 'Flye-contam']

quality_df = pd.merge(quality_df,polish_eval[['Mean_cov LR']],
                      left_index=True, right_index=True)

# Compare default model in all polishing steps
quality_df['Same_model'] = quality_df[['Polca-model', 'Poly-model', 'Medaka-model', 'Flye-model']].nunique(axis=1).apply(lambda x: 'Yes' if x == 1 else 'No')
quality_df['Dif_general_model'] = quality_df['Polca-compl_general'] - quality_df['Flye-compl_general']
quality_df['Dif_specific_model'] = quality_df['Polca-compl_specific'] - quality_df['Flye-compl_specific']

# Add Quality columns

ls_completeness = ['Polca-compl_general', 'Polca-compl_specific', 'Poly-compl_general','Poly-compl_specific',
                   'Medaka-compl_general', 'Medaka-compl_specific', 'Flye-compl_general','Flye-compl_specific']
ls_contamination = ['Polca-contam', 'Polca-contam', 'Poly-contam', 'Poly-contam',
                    'Medaka-contam', 'Medaka-contam', 'Flye-contam','Flye-contam']
ls_qual = ['Polca-qual-general', 'Polca-qual-specific', 'Poly-qual-general', 'Poly-qual-specific',
           'Medaka-qual-general', 'Medaka-qual-specific', 'Flye-qual-general', 'Flye-qual-specific']

for qual, compl, contam in zip(ls_qual, ls_completeness,ls_contamination):
    quality_df[qual] = 'low_quality' # default value to be updated
    quality_df.loc[(quality_df[compl] >= 50) & (quality_df[contam] <= 10), qual] = 'medium_quality'
    quality_df.loc[(quality_df[compl] >= 90) & (quality_df[contam] <= 5), qual] = 'high_quality'

# General model dataframe
quality_df_general = quality_df[quality_df['Polca-model']=='Gradient Boost (General Model)']
quality_df_general = quality_df_general.drop(columns=[col for col in quality_df_general.columns if 'specific' in col])
quality_df_general['Qual-evol'] = (
    quality_df_general['Flye-qual-general'] + '_' + quality_df_general['Polca-qual-general'])

quality_df_general['Qual-evol'] = quality_df_general['Qual-evol'].str.replace('high_quality_high_quality','High-to-High')
quality_df_general['Qual-evol'] = quality_df_general['Qual-evol'].str.replace('medium_quality_medium_quality','Medium-to-Medium')
quality_df_general['Qual-evol'] = quality_df_general['Qual-evol'].str.replace('medium_quality_high_quality','Medium-to-High')
quality_df_general['Qual-evol'] = quality_df_general['Qual-evol'].str.replace('low_quality_medium_quality','Low-to-Medium')
quality_df_general['Qual-evol'] = quality_df_general['Qual-evol'].str.replace('high_quality_medium_quality','High-to-Medium')

quality_df_general.columns = ['MIMAG-table_compl', 'Polca-compl', 'Polca-model', 'Polca-contam',
                               'Poly-compl', 'Poly-model', 'Poly-contam', 'Medaka-compl', 'Medaka-model',
                               'Medaka-contam', 'Flye-compl', 'Flye-model', 'Flye-contam', 'Mean_cov LR',
                               'Same_model', 'Diff-completeness', 'Polca-qual', 'Poly-qual', 'Medaka-qual',
                               'Flye-qual', 'Qual-evol']

quality_df_general = quality_df_general.reset_index()

# Specific model dataframe
quality_df_specific = quality_df[quality_df['Polca-model']=='Neural Network (Specific Model)']
quality_df_specific = quality_df_specific.drop(columns=[col for col in quality_df_specific.columns if 'general' in col])
quality_df_specific['Qual-evol'] = (
    quality_df_specific['Flye-qual-specific'] + '_' + quality_df_specific['Polca-qual-specific'])

quality_df_specific['Qual-evol'] = quality_df_specific['Qual-evol'].str.replace('high_quality_high_quality','High-to-High')
quality_df_specific['Qual-evol'] = quality_df_specific['Qual-evol'].str.replace('medium_quality_medium_quality','Medium-to-Medium')
quality_df_specific['Qual-evol'] = quality_df_specific['Qual-evol'].str.replace('medium_quality_high_quality','Medium-to-High')
quality_df_specific['Qual-evol'] = quality_df_specific['Qual-evol'].str.replace('low_quality_medium_quality','Low-to-Medium')
quality_df_specific['Qual-evol'] = quality_df_specific['Qual-evol'].str.replace('high_quality_medium_quality','High-to-Medium')

quality_df_specific.columns = ['MIMAG-table_compl', 'Polca-compl', 'Polca-model', 'Polca-contam',
                               'Poly-compl', 'Poly-model', 'Poly-contam', 'Medaka-compl', 'Medaka-model',
                               'Medaka-contam', 'Flye-compl', 'Flye-model', 'Flye-contam', 'Mean_cov LR',
                               'Same_model', 'Diff-completeness', 'Polca-qual', 'Poly-qual', 'Medaka-qual',
                               'Flye-qual', 'Qual-evol']

quality_df_specific = quality_df_specific.reset_index()

# Concatenate back specific and general model dataframes

quality_df_final = pd.concat([quality_df_general,quality_df_specific],axis=0)
quality_df_final = quality_df_final.set_index('index')
quality_df_final.to_csv('intermediate-outputs/polishing_evaluation/results/Checkm2_final_polish_eval.csv')