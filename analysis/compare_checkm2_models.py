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

# Calculate completeness differences
def calculate_completeness_diff(quality_df):
    final_diff = []
    final_diff_Flye_Med = []
    final_diff_Flye_Poly = []
    polca_compl = []
    poly_compl = []
    medaka_compl = []
    flye_compl = []
    qual_evol = []

    for index, row in quality_df.iterrows():
        if 'Specific' in row['Polca-model']:
            diff = row['Polca-compl_specific'] - row['Flye-compl_specific']
            flye_med_diff = row['Medaka-compl_specific'] - row['Flye-compl_specific']
            flye_poly_diff = row['Poly-compl_specific'] - row['Flye-compl_specific']

            polca_compl.append(row['Polca-compl_specific'])
            poly_compl.append(row['Poly-compl_specific'])
            medaka_compl.append(row['Medaka-compl_specific'])
            flye_compl.append(row['Flye-compl_specific'])

            qual_combination = row['Flye-qual-specific'] + '_' + row['Polca-qual-specific']

        elif 'General' in row['Polca-model']:
            diff = row['Polca-compl_general'] - row['Flye-compl_general']
            flye_med_diff = row['Medaka-compl_general'] - row['Flye-compl_general']
            flye_poly_diff = row['Poly-compl_general'] - row['Flye-compl_general']

            polca_compl.append(row['Polca-compl_general'])
            poly_compl.append(row['Poly-compl_general'])
            medaka_compl.append(row['Medaka-compl_general'])
            flye_compl.append(row['Flye-compl_general'])

            qual_combination = row['Flye-qual-general'] + '_' + row['Polca-qual-general']

        final_diff.append(diff)
        final_diff_Flye_Med.append(flye_med_diff)
        final_diff_Flye_Poly.append(flye_poly_diff)
        qual_evol.append(qual_combination)

    # Assign new columns
    quality_df['Flye_completeness'] = flye_compl
    quality_df['Medaka_completeness'] = medaka_compl
    quality_df['Polypolish_completeness'] = poly_compl
    quality_df['Polca_completeness'] = polca_compl
    quality_df['Diff_Flye-Med'] = final_diff_Flye_Med
    quality_df['Diff_Flye-Poly'] = final_diff_Flye_Poly
    quality_df['Diff_Flye-Polca'] = final_diff
    quality_df['Qual-evol'] = qual_evol

    # Improve Quality Evolution Labels
    replacements = {
        'high_quality_high_quality': 'High-to-High',
        'medium_quality_medium_quality': 'Medium-to-Medium',
        'medium_quality_high_quality': 'Medium-to-High',
        'low_quality_medium_quality': 'Low-to-Medium',
        'high_quality_medium_quality': 'High-to-Medium'
    }
    quality_df['Qual-evol'].replace(replacements, inplace=True)

    return quality_df

quality_df_final = calculate_completeness_diff(quality_df)

# Diff in contamination value after polishing steps
quality_df_final['Diff_contam_Flye_Med'] = quality_df_final['Medaka-contam'] - quality_df_final['Flye-contam']
quality_df_final['Diff_contam_Flye_Poly'] = quality_df_final['Poly-contam'] - quality_df_final['Flye-contam']
quality_df_final['Diff_contam_Flye_Polca'] = quality_df_final['Polca-contam'] - quality_df_final['Flye-contam']

# Re-order and clean-up unused columns
desired_order = ['MIMAG-table_compl', 'Mean_cov LR', 'Final_model', 'Flye_completeness', 'Medaka_completeness',
                 'Polypolish_completeness', 'Polca_completeness', 'Flye-contam', 'Medaka-contam',
                 'Poly-contam', 'Polca-contam', 'Diff_Flye-Med', 'Diff_Flye-Poly',
                 'Diff_Flye-Polca', 'Diff_contam_Flye_Med', 'Diff_contam_Flye_Poly',
                 'Diff_contam_Flye_Polca', 'Qual-evol']
quality_df_final = quality_df_final[desired_order]

# We chose the Checkm2 'default' model for the most polished MAG (Polca) as the reference
quality_df_final.to_csv('intermediate-outputs/polishing_evaluation/results/Checkm2_final_polish_eval.csv')
