import os
import pandas as pd

os.chdir('/data/Projects/ShanghaiDogs/')
ngs_raw_counts = pd.read_csv('intermediate-outputs/00_quality_control/ngs_raw_counts.csv',\
                             index_col=0)

ngs_raw_counts['bps']=ngs_raw_counts['Count']*150
ngs_raw_counts['sample']=ngs_raw_counts['File'].str.replace('.fq1.gz','')
ngs_raw_counts['sample']=ngs_raw_counts['sample'].str.replace('.fq2.gz','')

sample_summary = ngs_raw_counts.groupby('sample').sum()
sample_summary.drop('File',axis=1,inplace=True)
sample_summary.drop('Count',axis=1,inplace=True)
sample_summary = sample_summary.reset_index()

for index, row in sample_summary.iterrows():
    filename = f"intermediate-outputs/singlem_profiling/{row['sample'][:4]}/{row['sample']}_size.tsv"
    print(filename)
    with open(filename, 'w') as file:
        file.write("sample\tnum_bases\n")
        file.write(f"{row['sample']}\t{row['bps']}\n")