# Required libraries to install:
# pip install pandas matplotlib seaborn tabulate openpyxl

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from tabulate import tabulate

# Updated file paths
mags_report = "../data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv"
args_mag = "../intermediate-outputs/06_ARG/MAGs-ARGs_ALL_filt.txt"

df1 = pd.read_csv(mags_report)
df2 = pd.read_csv(args_mag)

# Merge the DataFrames based on the "Bin ID" column
merged_df = pd.merge(df1, df2, on="Bin ID", how="outer")

# Step 1: Extract species from the "Classification" column
def extract_species(classification):
    if pd.isna(classification):
        return "Unknown"
    parts = classification.split(';')
    for part in parts:
        if part.startswith('g__'):  # Identify genus
            genus = part[3:]
        if part.startswith('s__'):  # Identify species
            species = part[3:]
    if species:
        return species
    elif genus:  # If 's__' consider it as genus + "novel sp"
        return f"{genus} novel_sp"
    else:
        return "Unknown"  # If neither genus nor species are found

merged_df['Species'] = merged_df['Classification'].apply(extract_species)

# Step 2: Count ARGs per Bin ID in the merged data
arg_counts = merged_df.groupby('Bin ID').size().reset_index(name='ARG_Count')

# Step 3: Merge ARG counts back, keeping unique Bin ID-Species pairs
bin_species = merged_df[['Bin ID', 'Species']].drop_duplicates()
merged_with_counts = pd.merge(bin_species, arg_counts, on="Bin ID", how="left")

# Step 4: Handle cases where ARG_Count is NaN
merged_with_counts['ARG_Count'] = merged_with_counts['ARG_Count'].fillna(0)

# Step 5-8: Compute and display stats for all ARG counts (instead of saving to Excel)
overall_stats = pd.DataFrame({
    'Species': ['Overall'],
    'Mean': [merged_with_counts['ARG_Count'].mean()],
    'Median': [merged_with_counts['ARG_Count'].median()],
    'Min': [merged_with_counts['ARG_Count'].min()],
    'Max': [merged_with_counts['ARG_Count'].max()]
})

species_stats = merged_with_counts.groupby('Species')['ARG_Count'].agg(['mean', 'median', 'min', 'max']).reset_index()
species_stats.columns = ['Species', 'Mean', 'Median', 'Min', 'Max']

stats_df = pd.concat([overall_stats, species_stats], ignore_index=True)

# Format numbers
stats_df['Mean'] = stats_df['Mean'].round(6)
stats_df['Median'] = stats_df['Median'].round(3)
stats_df['Min'] = stats_df['Min'].round(2)
stats_df['Max'] = stats_df['Max'].round(2)

# Display stats in a table using tabulate
print("\nARG Statistics per Species:")
print("===========================")
print(tabulate(stats_df, headers='keys', tablefmt='pretty', showindex=False))

# Step 10: Prepare data for plotting
top_species = species_stats.sort_values(by='Mean', ascending=False).head(14)['Species']
df_top = merged_with_counts[merged_with_counts['Species'].isin(top_species)]

# Directory to save plots
plot_dir = "figures"
os.makedirs(plot_dir, exist_ok=True)

# Define a palette with exactly 14 colors
extended_palette = sns.color_palette("Dark2", n_colors=14)

# Step 11: Save Box Plots
sizes = [(4, 3), (5, 3)]
plot_filenames = ["boxplot_small.png", "boxplot_large.png"]

for size, filename in zip(sizes, plot_filenames):
    plt.figure(figsize=size)
    sns.boxplot(
        x='Species', 
        y='ARG_Count', 
        data=df_top, 
        hue='Species', 
        palette=extended_palette,
        dodge=False
    )
    
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.yticks(range(0, 61, 10), fontsize=8)
    plt.title('Boxplot of ARG Count Distribution', fontsize=10)
    plt.xlabel("Species", fontsize=9)
    plt.ylabel("Total ARG Count", fontsize=9)
    plt.legend([], [], frameon=False)
    
    plt.tight_layout()
    save_path = os.path.join(plot_dir, filename)
    plt.savefig(save_path, dpi=300)
    print(f"Saved: {save_path}")
    plt.close()

# Step 12: Save Violin Plot
violin_plot_path = os.path.join(plot_dir, "violin_plot.png")

plt.figure(figsize=(6, 4))
sns.violinplot(
    x='Species', 
    y='ARG_Count', 
    data=df_top, 
    hue='Species', 
    palette=extended_palette,
    inner="box",
    legend=False
)

plt.xticks(rotation=45, ha="right", fontsize=8)
plt.yticks(range(0, 61, 10), fontsize=8)
plt.title('Violin Plot of ARG Count Distribution', fontsize=10)
plt.xlabel("Species", fontsize=9)
plt.ylabel("Total ARG Count", fontsize=9)

plt.tight_layout()
plt.savefig(violin_plot_path, dpi=300)
print(f"Saved: {violin_plot_path}")
plt.close()

# Step 13: Find the top 10 ARGs by frequency
top_args = merged_df['Best_Hit_ARO'].value_counts().reset_index()
top_args.columns = ['Best_Hit_ARO', 'Count']  # Rename columns to match your example
top_10_args = top_args.head(10)  # Select top 10

# Display the results instead of saving
print("\nTop 10 ARGs by frequency:")
print(tabulate(top_10_args, headers='keys', tablefmt='pretty', showindex=False))

# Step 14: Find the top 5 most frequent ARGs per species
arg_counts_per_species = merged_df.groupby(['Species', 'Best_Hit_ARO']).size().reset_index(name='Counts')

# Sort by Species and Counts, then group by Species to get top 5 ARGs
top_args_per_species = arg_counts_per_species.sort_values(by=['Species', 'Counts'], ascending=[True, False])
top_5_args_per_species = top_args_per_species.groupby('Species').head(5).reset_index(drop=True)

# Rename columns to match your ideal output
top_5_args_per_species.columns = ['Species', 'Most Frequent ARGs', 'Counts']

# Save to Excel
top_args_species_output_path = "/work/microbiome/shanghai_dogs/output/Top_5_ARGs_per_Species.xlsx"
os.makedirs(os.path.dirname(top_args_species_output_path), exist_ok=True)  # Create output directory if it doesn't exist
top_5_args_per_species.to_excel(top_args_species_output_path, index=False)
print(f"Top 5 ARGs per Species saved to: {top_args_species_output_path}")
