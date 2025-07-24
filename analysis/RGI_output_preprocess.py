import sys
import pandas as pd
import os

os.chdir('/data/Projects/ShanghaiDogs')

# Import files

# MIMAG report
MIMAG_report = pd.read_csv("data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv", \
                  delimiter=',',header=0)

# concatenated RGI results on MAGs
arg_MAGs = pd.read_csv("intermediate-outputs/06_ARG/00_RGI_CARD/ALL_samples_strict_matched.txt", \
                       delimiter='\t',header=None,index_col=None)

# concatenated RGI results on contigs
arg_contigs = pd.read_csv("intermediate-outputs/06_ARG/01_RGI_CARD_contigs/contigs-ARGs_ALL_samples.txt", \
                          delimiter='\t',header=None,index_col=None)


# Reformat RGI output tables

## 1) For arg_MAGs
arg_MAGs.columns = ["ORF_ID", "Contig", "Start", "Stop", "Orientation", "Cut_Off", "Pass_Bitscore", "Best_Hit_Bitscore",
                    "Best_Hit_ARO", "Best_Identities", "ARO", "Model_type", "SNPs_in_Best_Hit_ARO", "Other_SNPs",
                    "Drug Class", "Resistance Mechanism", "AMR Gene Family", "Predicted_DNA", "Predicted_Protein",
                    "CARD_Protein_Sequence", "Percentage Length of Reference Sequence", "ID", "Model_ID", "Nudged", "Note",
                    "Hit_Start", "Hit_End", "Antibiotic"]

arg_MAGs['sample_id'] = arg_MAGs['ORF_ID'].str.extract(r'(\bD[0-9]+\b)')
arg_MAGs['bin_id'] = arg_MAGs['ORF_ID'].str.extract(r'(out_[^\/]+_SemiBin_[^\.]+)')
arg_MAGs['bin_id'] = arg_MAGs['bin_id'].str.replace('out_','')
arg_MAGs['bin_id'] = arg_MAGs['bin_id'] + '_' + arg_MAGs['sample_id']
arg_MAGs['partial']=arg_MAGs['ORF_ID'].str.extract(r'partial=(\d{2})')

arg_MAGs_rf=pd.merge(arg_MAGs,MIMAG_report[['Original ID','Quality','Bin ID']],left_on='bin_id',right_on='Original ID')

# Stricter filtering
arg_MAGs_filt = arg_MAGs_rf[arg_MAGs_rf['Best_Identities']>=90]
arg_MAGs_filt = arg_MAGs_filt[arg_MAGs_filt['Percentage Length of Reference Sequence']>=90]

arg_MAGs_filt.to_csv('intermediate-outputs/06_ARG/SHD1_ARGs_MAGs.csv.gz', index=False, compression='gzip')

## 2) For arg_contigs
arg_contigs.columns = ["ORF_ID", "Contig", "Start", "Stop", "Orientation", "Cut_Off", "Pass_Bitscore", "Best_Hit_Bitscore",
                    "Best_Hit_ARO", "Best_Identities", "ARO", "Model_type", "SNPs_in_Best_Hit_ARO", "Other_SNPs",
                    "Drug Class", "Resistance Mechanism", "AMR Gene Family", "Predicted_DNA", "Predicted_Protein",
                    "CARD_Protein_Sequence", "Percentage Length of Reference Sequence", "ID", "Model_ID", "Nudged", "Note",
                    "Hit_Start", "Hit_End", "Antibiotic", "sample_id"]

arg_contigs['partial']=arg_contigs['ORF_ID'].str.extract(r'partial=(\d{2})')

# Format is weird, so I extract the information from the ORF_id column
arg_contigs['Contig'] = arg_contigs['ORF_ID'].str.extract(r'(contig_[^ ]*\d+)')
arg_contigs['Start'] = arg_contigs['ORF_ID'].str.extract(r'#\s*(\d+)')
arg_contigs['Stop'] = arg_contigs['ORF_ID'].str.extract(r'#\s*\d+\s*#\s*(\d+)')
arg_contigs['Orientation'] = arg_contigs['ORF_ID'].str.extract(r'#\s*\d+\s*#\s*\d+\s*#\s*([+-])')
arg_contigs['Orientation'] = arg_contigs['Orientation'].fillna('+')

# Stricter filtering
arg_contigs_filt = arg_contigs[arg_contigs['Best_Identities']>=90]
arg_contigs_filt = arg_contigs_filt[arg_contigs_filt['Percentage Length of Reference Sequence']>=90]

arg_contigs_filt.to_csv('intermediate-outputs/06_ARG/SHD1_ARGs_contigs.csv.gz', index=False, compression='gzip')
