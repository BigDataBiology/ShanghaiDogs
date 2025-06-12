import pandas as pd

# paths
mimag_report_path = 'shanghai_dogs/data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv'
args_filt_path = 'microbiome/shanghai_dogs/intermediate-outputs/06_ARG/MAGs-ARGs_ALL_filt.txt'
mimag_df = pd.read_csv(mimag_report_path)
args_df = pd.read_csv(args_filt_path)

# extract species
mimag_df['Species'] = mimag_df['Classification'].str.extract(r's__([^;]+)')

# table with all Bin IDs and species
base_df = mimag_df[['Bin ID', 'Species']].copy()

# count ARGs per Bin
arg_counts = args_df.pivot_table(
    index='Bin ID',
    columns='Best_Hit_ARO',
    aggfunc='size',
    fill_value=0
).reset_index()

# merge with ARG counts 
full_df = pd.merge(base_df, arg_counts, on='Bin ID', how='left')

# fbins with no ARGs with 0
full_df = full_df.fillna(0)

#  columns to int
for col in full_df.columns[2:]:
    full_df[col] = full_df[col].astype(int)

# total ARGs column
full_df['Total_ARGs'] = full_df.iloc[:, 2:].sum(axis=1)

# save
full_df.to_csv('ARGs_by_Bin_Summary.csv', index=False)

print("File saved as 'ARGs_by_Bin_Summary.csv'")
