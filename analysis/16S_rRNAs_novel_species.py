# -*- coding: utf-8 -*-

"""
Created on Mon May 13 11:42:37 2024
@author: Anna Cusco
"""

import os
import pandas as pd

os.chdir('/data/Projects/ShanghaiDogs/')
SHD_table = pd.read_csv("/data/Projects/ShanghaiDogs/data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv")

# Keep novel species/taxonomy only
SHD_novel = SHD_table[(SHD_table['Representative'] == 'Yes') & (SHD_table['Classification'].str.endswith('_'))]

# Create a list of the 'novel' species
SHD_novel_ls = SHD_novel['Bin ID'].to_list() # new nomenclature
SHD_novel_ls = [('_'.join(id.split('_')[:-1])) for id in SHD_novel['Original ID']] # old nomenclature
SHD_novel_ls_sampleid = SHD_novel['Sample'].to_list()
files_list = [f"{x}/{x}_{y}_ribosomal.fa" for x, y in zip(SHD_novel_ls_sampleid, SHD_novel_ls)]

# Create a list with the full path to rRNA fasta files
rRNA_path = '/data/Projects/ShanghaiDogs/intermediate-outputs/07_ribosomal_genes/barrnap_fasta/'

all_16S_lines = []

for file_name in files_list:
    file_path = os.path.join(rRNA_path, file_name)
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i in range(len(lines)):
            if '16S' in lines[i]:
                all_16S_lines.append('>' + file_name + '--' + lines[i].strip())
                if i+1 < len(lines):
                    all_16S_lines.append(lines[i+1].strip()) # fasta sequence

# Define the file path where you want to store the data
output_file_path = '/data/Projects/ShanghaiDogs/intermediate-outputs/07_ribosomal_genes/16SrRNA_novel_sp.fa'

# Open the file in write mode
with open(output_file_path, 'w') as file:
    for line in all_16S_lines:
        file.write(line + '\n')