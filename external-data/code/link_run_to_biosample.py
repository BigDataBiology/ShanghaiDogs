import pandas as pd
import os

os.chdir('/data/Projects/ShanghaiDogs/external-data/')
SRA_metadata = pd.read_table('data/dog_microbiome_archive_otu_tables/SRR_Acc_List_Metadata.txt', sep='|',\
                             skiprows=[1],index_col=0)

columns = list(SRA_metadata.columns)
run_to_biosample = SRA_metadata[[' biosample      ']]
run_to_biosample = run_to_biosample.reset_index()
run_to_biosample.columns = ['run','biosample']
run_to_biosample['run'] = run_to_biosample['run'].str.replace(' ', '', regex=True)
run_to_biosample['biosample'] = run_to_biosample['biosample'].str.replace(' ', '', regex=True)

# Create a dictionary to store lists of run IDs for each biosample
biosample_run_lists = {}

for run, biosample in run_to_biosample.values:
    if biosample not in biosample_run_lists:
        biosample_run_lists[biosample] = []
    biosample_run_lists[biosample].append(run)

# Write the lists of run numbers to separate biosample.txt files
for biosample, run_list in biosample_run_lists.items():
    with open(f'data/dog_microbiome_archive_otu_tables/run_to_biosample/{biosample}_runs_list.txt', 'w') as f:
        for run_number in run_list:
            full_path=f"../renew_outputs/{run_number[:5]}/{run_number}.json"
            f.write(f"{full_path}\n")