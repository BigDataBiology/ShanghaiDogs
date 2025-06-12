import pandas as pd

# paths
mimag_report_path = '/work/microbiome/shanghai_dogs/data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv'
args_filt_path = '/work/microbiome/shanghai_dogs/intermediate-outputs/06_ARG/MAGs-ARGs_ALL_filt.txt'

# input files
mimag_df = pd.read_csv(mimag_report_path)
args_df = pd.read_csv(args_filt_path)

# Merge the dataframes using 'Bin ID'
merged_df = pd.merge(mimag_df, args_df, on='Bin ID', how='left')

# pivot table to count occurrences of each ARG per Bin ID
arg_counts = merged_df.pivot_table(
    index='Bin ID',
    columns='Best_Hit_ARO',
    aggfunc='size',
    fill_value=0
)

# make 'Bin ID' a column again
arg_counts = arg_counts.reset_index()

arg_counts['Total_ARGs'] = arg_counts.iloc[:, 1:].sum(axis=1)

# Save table
arg_counts.to_csv('ARGs_by_Bin_Summary.csv', index=False)

print("table saved as 'ARGs_by_Bin_Summary.csv'")

