# -*- coding: utf-8 -*-

"""
Created on Tue Feb 20 13:43:05 2024
@author: Anna Cusco
"""

import sys
import pandas as pd
import os
import gzip

sample_id = sys.argv[1]
#sample_id='D000'

# Import checkm2, gtdbtk and contig count tables
checkm = "/data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool/"+sample_id+"/checkm2/quality_report.tsv"
gtdbtk = "/data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool/"+sample_id+"/gtdbtk/gtdbtk.bac120.summary.tsv"
tigs = "/data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool/"+sample_id+"/tigs_num_"+sample_id+".txt"

with open(checkm, 'r') as file:
  checkm = pd.read_csv(file, delimiter='\t',header=0,index_col=0)
  checkm_sub = checkm[['Completeness','Contamination','Genome_Size']]
  checkm_sub.index=checkm_sub.index.str.replace('.fa.gz','')
  suffix='_'+sample_id
  checkm_sub.index=checkm_sub.index.str.replace(suffix,'')

with open(gtdbtk, 'r') as file:
  gtdbtk = pd.read_csv(file, delimiter='\t',header=0,index_col=0)
  gtdbtk_sub = gtdbtk[['classification','fastani_reference']]
  suffix='_'+sample_id
  gtdbtk_sub.index=gtdbtk_sub.index.str.replace(suffix,'')

with open(tigs, 'r') as file:
  tigs = pd.read_csv(file, delimiter=',',header=0,index_col=0)
  tigs.index = tigs.index.str.replace('.fa', '')
  suffix = '_' + sample_id
  tigs.index=tigs.index.str.replace(suffix,'')

# Create 'Quality' column
qual_report = pd.merge(checkm_sub,gtdbtk_sub,left_index=True,right_index=True)
qual_report = pd.merge(qual_report,tigs,left_index=True,right_index=True)

qual_report['Quality'] = 'low_quality'  # Default value
qual_report.loc[(qual_report['Completeness'] >= 50) & (qual_report['Contamination'] < 10), 'Quality'] = 'medium_quality'
qual_report.loc[(qual_report['Completeness'] >= 90) & (qual_report['Contamination'] < 5), 'Quality'] = 'high_quality'
qual_report['sample_id']=sample_id

qual_report['complete_ID']=qual_report.index+'_'+qual_report['sample_id']
qual_report=qual_report.reset_index()
qual_report.columns = ['Bin ID','Completeness','Contamination','Genome_Size','Classification','fastani_reference',\
                       'Nr contigs', 'Quality', 'Sample', 'Complete ID']

# Add ribosomal genes
rib = "/data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool/"+sample_id+"/barrnap_out/"+sample_id+"_ALL.txt"

with open(rib, 'r') as file:
  rib = pd.read_csv(file, delimiter='\t',header=None)
  rib.columns = ["file_name","seqname","source","feature","start","end","score","strand","frame","attribute"]
  name=sample_id+"_summary.txt"
  rib=rib[rib['file_name'].str.contains(name) == False]
  rib['sample_id'] = sample_id
  rib['bin_id'] = rib['file_name'].str.replace(sample_id+'_','')
  rib['bin_id'] = rib['bin_id'].str.replace('_barrnap.txt', '')
  rib['Ribosomal gene'] = rib['attribute'].str.replace(r'^.*product=', '', regex=True)
  rib['Ribosomal gene']= rib['Ribosomal gene'].str.replace(r'\(partial\).*$','partial',regex=True)

pivot_rib = rib.pivot_table(index=['bin_id'], columns='Ribosomal gene',aggfunc='size', fill_value=0)

if '16S ribosomal RNA partial' not in pivot_rib.columns:
  pivot_rib['16S ribosomal RNA partial'] = 0

if '23S ribosomal RNA partial' not in pivot_rib.columns:
  pivot_rib['23S ribosomal RNA partial'] = 0

if '5S ribosomal RNA partial' not in pivot_rib.columns:
  pivot_rib['5S ribosomal RNA partial'] = 0

pivot_rib = pivot_rib[['16S ribosomal RNA','23S ribosomal RNA','5S ribosomal RNA',\
                       '16S ribosomal RNA partial','23S ribosomal RNA partial','5S ribosomal RNA partial']]
pivot_rib.columns = ['16S rRNA','23S rRNA','5S rRNA','16S partial','23S partial','5S partial']

qual_report_rib = pd.merge(qual_report, pivot_rib,left_on='Bin ID',right_index=True,how='outer')

dest_folder = "/data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool/"+sample_id+"/"+sample_id+"_qual_report.csv"
qual_report_rib.to_csv(dest_folder, index=False)
print(f"DataFrame has been saved to '{dest_folder}'")

# Add ARGs
arg = "/data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool/"+sample_id+"/RGI_CARD/"+sample_id+"_ARG_strict_matches.txt.gz"
arg = pd.read_csv(arg, compression='gzip', delimiter='\t', header=None)

arg.columns = ["ORF_ID", "Contig", "Start", "Stop", "Orientation", "Cut_Off", "Pass_Bitscore", "Best_Hit_Bitscore",
               "Best_Hit_ARO", "Best_Identities", "ARO", "Model_type", "SNPs_in_Best_Hit_ARO", "Other_SNPs",
               "Drug Class", "Resistance Mechanism", "AMR Gene Family", "Predicted_DNA", "Predicted_Protein",
               "CARD_Protein_Sequence", "Percentage Length of Reference Sequence", "ID", "Model_ID", "Nudged", "Note",
               "Hit_Start", "Hit_End", "Antibiotic"]
arg['sample_id'] = arg['ORF_ID'].str.extract(r'(\bD[0-9]+\b)')
arg['bin_id'] = arg['ORF_ID'].str.extract(r'(out_[^\/]+_SemiBin_[^\.]+)')
arg['bin_id'] = arg['bin_id'].str.replace('out_','')
arg['partial']=arg['ORF_ID'].str.extract(r'partial=(\d{2})')
arg=arg[arg['Best_Identities']>=90]
arg=arg[arg['Percentage Length of Reference Sequence']>=90]

pivot_arg = arg.pivot_table(index=['bin_id'], columns='Best_Hit_ARO',aggfunc='size', fill_value=0)
pivot_arg['Total_ARGs']=pivot_arg.sum(axis=1)

qual_report_arg = pd.merge(qual_report_rib, pivot_arg['Total_ARGs'],left_on='Bin ID',right_index=True,how='outer')
qual_report_arg=qual_report_arg.fillna(0)

dest_folder = "/data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool/"+sample_id+"/"+sample_id+"_qual_report_ARGs.csv"
qual_report_arg.to_csv(dest_folder, index=False)
print(f"DataFrame has been saved to '{dest_folder}'")

sys.exit()
