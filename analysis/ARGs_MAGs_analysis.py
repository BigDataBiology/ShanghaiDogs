import sys
import pandas as pd
import os

# Import files
MIMAG_report = pd.read_csv("/data/Projects/ShanghaiDogs/data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv", \
                  delimiter=',',header=0)

arg = pd.read_csv("/data/Projects/ShanghaiDogs/intermediate-outputs/06_ARG/00_RGI_CARD/ALL_samples_strict_matched.txt", \
                  delimiter='\t',header=None,index_col=None)

# Reformat RGI output table
arg.columns = ["ORF_ID", "Contig", "Start", "Stop", "Orientation", "Cut_Off", "Pass_Bitscore", "Best_Hit_Bitscore",
              "Best_Hit_ARO", "Best_Identities", "ARO", "Model_type", "SNPs_in_Best_Hit_ARO", "Other_SNPs",
              "Drug Class", "Resistance Mechanism", "AMR Gene Family", "Predicted_DNA", "Predicted_Protein",
              "CARD_Protein_Sequence", "Percentage Length of Reference Sequence", "ID", "Model_ID", "Nudged", "Note",
              "Hit_Start", "Hit_End", "Antibiotic"]

arg['sample_id'] = arg['ORF_ID'].str.extract(r'(\bD[0-9]+\b)')
arg['bin_id'] = arg['ORF_ID'].str.extract(r'(out_[^\/]+_SemiBin_[^\.]+)')
arg['bin_id'] = arg['bin_id'].str.replace('out_','')
arg['bin_id'] = arg['bin_id'] + '_' + arg['sample_id']
arg['partial']=arg['ORF_ID'].str.extract(r'partial=(\d{2})')

arg_MAGs=pd.merge(arg,MIMAG_report[['Original ID','Quality','Bin ID']],left_on='bin_id',right_on='Original ID')

# Stricter filtering
arg_MAGs = arg_MAGs[arg_MAGs['Best_Identities']>=90]
arg_MAGs = arg_MAGs[arg_MAGs['Percentage Length of Reference Sequence']>=90]

