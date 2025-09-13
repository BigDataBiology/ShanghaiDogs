# -*- coding: utf-8 -*-

"""
Created on Thu Apr 11 15:24:17 2024
@author: Anna Cusco
"""

import os
import glob
import pandas as pd

# Import checkm2, contig count

checkm = pd.read_csv('external-data/data/NCBI_genomes_ref/checkm2/quality_report.tsv',sep='\t')
checkm_sub = checkm[['Name','Completeness','Contamination','Genome_Size']]
checkm_sub['Quality'] = 'low-quality'  # Default value
checkm_sub.loc[(checkm_sub['Completeness'] >= 50) & (checkm_sub['Contamination'] < 10), 'Quality'] = 'medium-quality'
checkm_sub.loc[(checkm_sub['Completeness'] >= 90) & (checkm_sub['Contamination'] < 5), 'Quality'] = 'high-quality'

tigs = pd.read_csv('external-data/data/NCBI_genomes_ref/contigs_count.txt',sep=',')
tigs['Filename']=tigs['Filename'].str.replace('.fna.gz','')

qual_report = pd.merge(checkm_sub,tigs,left_on='Name',right_on='Filename',how='outer')
qual_report['Name'] = qual_report['Name'].apply(lambda x: '_'.join(x.split('_')[:2]))

# Add ribosomal genes
rRNA_path = 'external-data/data/NCBI_genomes_ref/barrnap/out/'
all_rRNAs = []

for filename in glob.glob(os.path.join(rRNA_path, '*.txt')):
  with open(filename, 'r') as file:
    content = file.read().strip()
    if content == "##gff-version 3":
      print("No ribosomal genes were detected using barrnap")
    else:
      print("Ribosomal genes detected... counting them")
      file.seek(0)
      rRNA = pd.read_csv(file, sep='\t', skiprows=[0],header=None)
      rRNA['Name']=filename
      rRNA['Name']=rRNA['Name'].str.replace(rRNA_path,'')
      rRNA['Name'] = rRNA['Name'].str.replace('_barrnap.txt', '')
      rRNA.columns = ["seqname","source","feature","start","end","score","strand","frame","attribute","Name"]
      rRNA['Ribosomal gene'] = rRNA['attribute'].str.replace(r'^.*product=', '', regex=True)
      rRNA['Ribosomal gene'] = rRNA['Ribosomal gene'].str.replace(r'\(partial\).*$', 'partial', regex=True)
      all_rRNAs.append(rRNA)

all_rRNAs_df = pd.concat(all_rRNAs, ignore_index=True)
pivot_rRNA=all_rRNAs_df.pivot_table(index=['Name'],columns='Ribosomal gene',aggfunc='size',fill_value=0)
pivot_rRNA = pivot_rRNA.reset_index()

# Total includes all detected, complete + partial
pivot_rRNA = pivot_rRNA[['Name','16S ribosomal RNA','23S ribosomal RNA','5S ribosomal RNA',\
                         '16S ribosomal RNA partial','23S ribosomal RNA partial','5S ribosomal RNA partial']]
pivot_rRNA.columns = ['Name','16S rRNA','23S rRNA','5S rRNA','16S partial','23S partial','5S partial']

qual_report = pd.merge(qual_report,pivot_rRNA,left_on='Name',right_on='Name',how='outer')
qual_report = qual_report.fillna(0)

# Add tRNAs

tRNA_path = 'external-data/data/NCBI_genomes_ref/tRNAs/'
tRNA_count = pd.DataFrame(columns=['Name', 'Unique tRNAs', 'Total tRNAs'])

all_tRNAs = []
for filename in glob.glob(os.path.join(tRNA_path, '*.out')):
  print(filename)
  tRNA_out = pd.read_csv(filename, skiprows=3, delimiter='\t',header=None)
  tRNA_out.columns = ("contig_id","tRNA","Begin","End","Type","Codon","Begin","End","Score","Note")
  tRNA_out['Name'] = filename
  tRNA_out['Name'] = tRNA_out['Name'].str.replace(tRNA_path, '')
  tRNA_out['Name'] = tRNA_out['Name'].str.replace('_trna.out','')
  aa_count = tRNA_out['Type'].value_counts() # Count each AA Type
  aa_count = aa_count.reset_index()
  unique_aa = len(aa_count['Type'])
  total_aa = aa_count['count'].sum()
  tRNA_count.loc[tRNA_out['Name'][0]] = [tRNA_out['Name'][0], unique_aa, total_aa]

qual_report_final = pd.merge(qual_report,tRNA_count,left_on='Name',right_on='Name',how='outer')

# Assess if the genomes on NCBI follow MIMAG criteria
qual_report_final['MIMAG']='No'

for index, row in qual_report_final.iterrows():
    if row['Quality'] == 'high-quality' and row['16S rRNA'] > 0 and row['23S rRNA'] > 0 and row['Unique tRNAs'] > 19:
            print(row['Quality'])
            qual_report_final.loc[index, 'MIMAG'] = 'Yes'

dest_folder = "external-data/data/NCBI_genomes_ref/NCBI_genomes_qual_MIMAG_report.csv"
qual_report_final.to_csv(dest_folder, index=False)
print(f"DataFrame has been saved to '{dest_folder}'")
