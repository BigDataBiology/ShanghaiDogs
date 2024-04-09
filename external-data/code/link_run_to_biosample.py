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

# Separate biosamples for number of runs
multiple_output_dir = 'data/dog_microbiome_archive_otu_tables/run_to_biosample/multiple_run'
single_output_dir = 'data/dog_microbiome_archive_otu_tables/run_to_biosample/single_run'

for biosample, run_list in biosample_run_lists.items():
    if len(run_list) == 1:
        output_file = os.path.join(single_output_dir, f'{biosample}_runs_list.txt')
    elif len(run_list) > 1:
        output_file = os.path.join(multiple_output_dir, f'{biosample}_runs_list.txt')
    with open(output_file, 'w') as ofile:
        for run_number in run_list:
            json_file = f"renew_outputs/{run_number[:5]}/{run_number}.json"
            if os.path.exists(f'data/dog_microbiome_archive_otu_tables/{json_file}'):
                ofile.write(f"../../{json_file}\n")
