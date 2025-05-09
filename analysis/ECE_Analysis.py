import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import gzip
import re

#constant
threshold = 0.8
IN2CM = 2.54

# Load metadata
metadata = pd.read_csv('/work/microbiome/shanghai_dogs/data/ShanghaiDogsMetadata/REP_canid_metadata.csv', sep=',', encoding='unicode_escape', index_col=0)
metadata['Sex'] = metadata['Sex'].str.lower()

# Load coverage and NCE info
MAGs_NCE_covered_frac = pd.read_csv('/work/microbiome/shanghai_dogs/intermediate-outputs/external_datasets_mappings/SHD_covered_fraction.tsv.gz',
                        sep='\t', index_col=0)
NCE_info = pd.read_csv('/work/microbiome/shanghai_dogs/data/ShanghaiDogsTables/SHD_NC_props.tsv.gz',
                       sep='\t', index_col=0)

# Load Contigs and Extrachromosomal elements data 
Contigs_file = "/work/microbiome/shanghai_dogs/intermediate-outputs/06_ARG/contigs-ARGs_ALL_filt.txt"
Extrachromosomal_elements_file = "/work/microbiome/shanghai_dogs/data/ShanghaiDogsTables/SHD_NC_props.tsv.gz"

df_arg = pd.read_csv(Contigs_file, sep=",")
df_non_chromo = pd.read_csv(Extrachromosomal_elements_file, sep="\t")

# Clean Contig column in df_arg
df_arg['Contig_clean'] = df_arg['Contig'].str.extract(r'(contig_\d+)')

# Merge dataframes
merged_df = pd.merge(
    df_non_chromo,
    df_arg[['Contig_clean', 'sample_id', 'Best_Hit_ARO', 'Cut_Off', 'Best_Identities']],
    left_on=['Contig', 'Sample'],
    right_on=['Contig_clean', 'sample_id'],
    how='left'
)

merged_df = merged_df.drop(columns=['Contig_clean'])

# Filter rows with valid Best_Hit_ARO and Strict/Perfect Cut_Off
filtered_df = merged_df[merged_df['Best_Hit_ARO'].notnull()]
filtered_df = filtered_df[filtered_df['Cut_Off'].isin(['Strict', 'Perfect'])]

# Create list of unique elements with AROs and categories
element_info_rows = filtered_df[['Element', 'Best_Hit_ARO', 'Category']].drop_duplicates()
element_info_list = list(element_info_rows.itertuples(index=False, name=None))

# Map abbreviations for metadata
abbreviations = {
    'Dog Pet': 'Current_Pet',
    'Dog Colony': 'Current_Colony',
    'Dog Shelter': 'Current_Shelter',
    'Dog Free_roaming': 'Current_Free'
}

# Calculate prevalence
prevalence_records = []
for element_id, aro, category in element_info_list:
    if element_id not in MAGs_NCE_covered_frac.index or element_id not in NCE_info.index:
        continue

    coverage = MAGs_NCE_covered_frac.loc[element_id]
    presence = coverage[coverage >= threshold]
    positive_samples = presence.index.str.split('_').str[0]

    metadata_env = metadata[['env_classification']].copy()
    metadata_env = metadata_env[metadata_env['env_classification'].isin(abbreviations.keys())]
    metadata_env['is_positive'] = metadata_env.index.isin(positive_samples)

    prevalence = metadata_env.groupby('env_classification')['is_positive'].mean() * 100

    record = {
        'Element': element_id,
        'Best_Hit_ARO': aro,
        'Category': category,
        'Current_Pet': prevalence.get('Dog Pet', np.nan),
        'Current_Colony': prevalence.get('Dog Colony', np.nan),
        'Current_Shelter': prevalence.get('Dog Shelter', np.nan),
        'Current_Free': prevalence.get('Dog Free_roaming', np.nan)
    }
    prevalence_records.append(record)

# Create dataframe
df = pd.DataFrame(prevalence_records)

# Extract size in kbp
df_non_chromo['Size'] = df_non_chromo['Working_header'].str.extract(r'Size_(\d+)').astype(float)
df_non_chromo['Size_kbp'] = df_non_chromo['Size'] / 1000
size_data = df_non_chromo[['Element', 'Size_kbp']].drop_duplicates()

# Merge size data
df = pd.merge(df, size_data, on='Element', how='left')

# Apply log transformation to sizes
df['Log_Size_kbp'] = np.log10(df['Size_kbp'])
min_log_size = df['Log_Size_kbp'].min()
max_log_size = df['Log_Size_kbp'].max()

# Normalize log sizes
df['Normalized_Log_Size'] = (df['Log_Size_kbp'] - min_log_size) / (max_log_size - min_log_size) * 0.8

# Load and process putative hosts
putative_hosts = pd.read_csv('/work/microbiome/shanghai_dogs/intermediate-outputs/external_datasets_mappings/dogs_putative_hosts.csv', sep=',')
putative_hosts['Species'] = putative_hosts['Putative host'].str.extract(r's__([^;]+)')
host_mapping = putative_hosts.groupby('Putative Plasmid')['Species'].first().to_dict()
df['Putative_Host'] = df['Element'].map(host_mapping)
df['Putative_Host'] = df['Putative_Host'].fillna('Unknown')

# Set up plot
n_rows = len(df)
fig_height = max(12 / IN2CM, n_rows * 0.7 / IN2CM)
fig = plt.figure(figsize=(12 / IN2CM, fig_height))
gs = gridspec.GridSpec(2, 6, height_ratios=[6, 0.5],
                       width_ratios=[0.5, 1.5, 1.5, 0.1, 0.9, 1.0], hspace=0.3, wspace=0.005)

# Define colormap
ylorbr_colors = cm.YlOrBr(np.linspace(0, 1, 256))
prevalence_cmap = LinearSegmentedColormap.from_list('YlOrBr_Custom', ylorbr_colors)

# Plot Category and Size (log scale)
ax0 = fig.add_subplot(gs[0, 0])
min_size_kbp = df['Size_kbp'].min()
max_size_kbp = df['Size_kbp'].max()
ax0.set_xlim(10, 500)
ax0.set_xscale('log')
ax0.set_ylim(0, n_rows)
ax0.set_xlabel('Size (kbp)', fontsize=9)
ax0.set_xticks([10, 50, 100])
ax0.set_xticklabels(['10', '50', '100'], fontsize=9)
ax0.set_yticks([])
ax0.spines['left'].set_visible(True)
ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax0.spines['bottom'].set_visible(True)

unique_categories = df['Category'].unique()
dark2_colors = sns.color_palette("Dark2", n_colors=len(unique_categories))
category_to_color = dict(zip(unique_categories, dark2_colors))
for i, (cat, size_kbp) in enumerate(zip(df['Category'], df['Size_kbp'])):
    ax0.barh(n_rows - i - 0.5, size_kbp, height=0.4, color=category_to_color[cat])
ax0.text(100, n_rows + 0.5, "Category\n& Size", fontsize=9, ha='center', rotation=45)

# Plot Element names
ax1 = fig.add_subplot(gs[0, 1])
ax1.set_xlim(0, 1)
ax1.set_ylim(0, n_rows)
ax1.axis('off')
for i, name in enumerate(df['Element']):
    ax1.text(0, n_rows - i - 0.5, name, va='center', fontsize=9)
ax1.text(0.5, n_rows + 0.5, "Element", fontsize=9, ha='center', rotation=45)

# Plot Heatmap
prevalence_columns = ['Current_Pet', 'Current_Colony', 'Current_Free', 'Current_Shelter']
full_titles = ['Pet', 'Colony', 'Free', 'Shelter']
heatmap_data = df[prevalence_columns].copy()
heatmap_data.columns = full_titles

ax2 = fig.add_subplot(gs[0, 2])
sns.heatmap(heatmap_data, ax=ax2, cmap=prevalence_cmap, vmin=0, vmax=50,
            cbar=False, xticklabels=True, yticklabels=False,
            linewidths=0, linecolor='none')

ax2.set_xticks(np.arange(len(full_titles)) + 0.5)
ax2.set_xticklabels(full_titles, rotation=45, ha='center', fontsize=9)
ax2.xaxis.tick_top()
ax2.set_xlabel('')

# Plot Colorbar
cb_ax = fig.add_subplot(gs[1, 2])
cb = plt.colorbar(ax2.collections[0], cax=cb_ax, orientation='horizontal')
cb.ax.tick_params(labelsize=9)

# Plot ARG labels
ax3 = fig.add_subplot(gs[0, 4])
ax3.set_xlim(0, 1)
ax3.set_ylim(0, n_rows)
ax3.axis('off')
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.spines['left'].set_visible(False)
perfect_aros = set(df_arg[df_arg['Cut_Off'] == 'Perfect']['Best_Hit_ARO'].unique())
strict_aros = set(df_arg[df_arg['Cut_Off'] == 'Strict']['Best_Hit_ARO'].unique())
for i, arg in enumerate(df['Best_Hit_ARO']):
    weight = 'bold' if arg in perfect_aros else ('normal' if arg in strict_aros else 'light')
    ax3.text(0, n_rows - i - 0.5, arg, va='center', fontsize=9)
ax3.text(0, n_rows + 0.5, "ARG", fontsize=9, ha='left', rotation=45)

# Plot Putative Host
ax4 = fig.add_subplot(gs[0, 5])
ax4.set_xlim(0, 1)
ax4.set_ylim(0, n_rows)
ax4.axis('off')
ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)
ax4.spines['bottom'].set_visible(False)
ax4.spines['left'].set_visible(False)
for i, host in enumerate(df['Putative_Host']):
    ax4.text(0, n_rows - i - 0.5, host, va='center', fontsize=9)
ax4.text(0, n_rows + 0.5, "Putative Host", fontsize=9, ha='left', rotation=45)

fig.tight_layout()
fig.subplots_adjust(top=0.95, bottom=0.10, hspace=0.3, wspace=0.005)

fig.savefig('prevalence_heatmap.png', bbox_inches='tight', dpi=300)

# Part 2: Circular Gene Plot 
faa_file = "/work/microbiome/shanghai_dogs/intermediate-outputs/Prodigal/D003/D003_proteins.faa.gz"
annotations_file = "/work/microbiome/shanghai_dogs/intermediate-outputs/eggNOG_annot_contigs/D003/D003.emapper.annotations"

# COG categories to consider
valid_cogs = set("DLISJLV")

# Step 1: Parse gene coordinates from .faa.gz
pattern = re.compile(r"#\s+(\d+)\s+#\s+(\d+)")
genes = []

with gzip.open(faa_file, "rt") as f:
    for line in f:
        if line.startswith(">contig_2645"):
            match = pattern.search(line)
            if match:
                start, end = map(int, match.groups())
                contig_id = line.split()[0][1:]  
                strand = int(line.split("#")[3].strip().split()[0])  
                genes.append({
                    "start": start,
                    "end": end,
                    "contig": contig_id,
                    "strand": strand
                })

# Step 2: Parse COG categories from annotations
cog_categories = {}

with open(annotations_file, "r") as f:
    for line in f:
        if line.startswith("#") or not line.startswith("contig_2645"):
            continue
        parts = line.strip().split("\t")
        if len(parts) > 6:
            contig_name = parts[0]
            cog_field = parts[6]
            filtered = "".join(c for c in cog_field if c in valid_cogs)
            cog_categories[contig_name] = filtered if filtered else "S"

# Step 3: Map COGs to genes
for gene in genes:
    gene["cog_category"] = cog_categories.get(gene["contig"], "S")

# Step 4: Plot
genomic_size = max(g["end"] for g in genes)
fig, ax = plt.subplots(figsize=(3, 3))  
fig_circle = plt.Circle((0, 0), 0.15, fill=False, color='black', linewidth=1.5)
ax.add_artist(fig_circle)

# Color map: Grey for S, Dark2 for others
unique_cogs = sorted(set(g["cog_category"] for g in genes))
color_map = {"S": [0.6, 0.6, 0.6]}  # Custom grey color (dark grey)
non_s_cogs = [c for c in unique_cogs if c != "S"]
color_map.update(zip(non_s_cogs, [c for i, c in enumerate(plt.cm.Dark2.colors) if i != 7][:len(non_s_cogs)]))

# Arcs and labels with strand indication (+/-)
for gene in genes:
    angle = ((gene["start"] + gene["end"]) / 2 / genomic_size) * 2 * np.pi
    arc_len = (gene["end"] - gene["start"]) / genomic_size * 2 * np.pi
    theta = np.linspace(angle - arc_len / 2, angle + arc_len / 2, 100)
    r = 0.15
    x, y = r * np.cos(theta), r * np.sin(theta)
    cog = gene["cog_category"]
    color = color_map.get(cog)

    # Reverse arc direction if strand is -1 (counterclockwise)
    if gene["strand"] == -1:
        x, y = x[::-1], y[::-1]

    ax.plot(x, y, linewidth=4, color=color)

    # Label with strand indication (+ for 1, - for -1)
    strand_symbol = "+" if gene["strand"] == 1 else "-"
    label = f"{cog}{strand_symbol}" if cog != "NA" else f"NA{strand_symbol}"
    ax.text(0.19 * np.cos(angle), 0.19 * np.sin(angle), label, ha='center', va='center', fontsize=9, color='black')

ax.set_aspect('equal')
ax.set_xlim(-0.3, 0.3)  
ax.set_ylim(-0.3, 0.3) 
ax.axis("off")

fig.tight_layout()
fig.savefig('circular_gene_plot.png')
