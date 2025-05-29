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
from Bio import SeqIO

# Constants
IN2CM = 2.54
threshold = 0.8

# Load data
metadata = pd.read_csv('data/ShanghaiDogsMetadata/REP_canid_metadata.csv', sep=',', encoding='unicode_escape', index_col=0)
MAGs_ECE_covered_frac = pd.read_csv('intermediate-outputs/external_datasets_mappings/SHD_covered_fraction.tsv.gz', sep='\t', index_col=0)
df_Extra_chromo = pd.read_csv('data/ShanghaiDogsTables/SHD1_EC_props.tsv.gz', sep='\t', index_col=0)
df_Extra_chromo = df_Extra_chromo.reset_index()  # Reset index to make 'Element' a column
df_arg = pd.read_csv('intermediate-outputs/06_ARG/contigs-ARGs_ALL_filt.txt', sep=",")

# Clean Contig column in df.arg
df_arg['Contig_clean'] = df_arg['Contig'].str.extract(r'(contig_\d+)')

# Merge dataframes
merged_df = pd.merge(
    df_Extra_chromo,
    df_arg[['Contig_clean', 'sample_id', 'Best_Hit_ARO', 'Cut_Off', 'Best_Identities']],
    left_on=['Contig', 'Sample'],
    right_on=['Contig_clean', 'sample_id'],
    how='left'
)
merged_df = merged_df.drop(columns=['Contig_clean'])

# Filter rows with valid Best_Hit_ARO and Strict/Perfect Cut_Off
filtered_df = merged_df[merged_df['Best_Hit_ARO'].notnull()]
filtered_df = filtered_df[filtered_df['Cut_Off'].isin(['Strict', 'Perfect'])]

# Create list of unique elements with AROs, categories, and Cut_Off
element_info_rows = filtered_df[['Element', 'Best_Hit_ARO', 'Category', 'Cut_Off']].drop_duplicates()
element_info_rows['Display_Element'] = element_info_rows['Element'].str.replace('NC', 'EC', regex=False)
element_info_list = list(element_info_rows[['Display_Element', 'Best_Hit_ARO', 'Category', 'Cut_Off']].itertuples(index=False, name=None))

# Calculate prevalence, splitting Pet into This_Study (Shanghai Dogs) and Current_Pet (other Pets)
prevalence_records = []
for display_element, aro, category, cutoff in element_info_list:
    element_id = display_element.replace('EC', 'NC')
    if element_id not in MAGs_ECE_covered_frac.index:
        continue
    coverage = MAGs_ECE_covered_frac.loc[element_id]
    presence = coverage[coverage >= threshold]
    positive_samples = presence.index.str.split('_').str[0]
    metadata_env = metadata[['env_classification', 'Study']].copy()
    metadata_env = metadata_env[metadata_env['env_classification'].isin(['Dog Pet', 'Dog Colony', 'Dog Shelter', 'Dog Free_roaming'])]
    metadata_env['is_positive'] = metadata_env.index.isin(positive_samples)
    
    # Split Dog Pet into Shanghai Dogs (This_study) and other Pets
    pet_samples = metadata_env[metadata_env['env_classification'] == 'Dog Pet']
    shanghai_pet = pet_samples[pet_samples['Study'] == 'This_study']
    other_pet = pet_samples[pet_samples['Study'] != 'This_study']
    
    # Calculate prevalence for each group
    prevalence = metadata_env.groupby('env_classification')['is_positive'].mean() * 100
    shanghai_pet_prev = shanghai_pet['is_positive'].mean() * 100 if not shanghai_pet.empty else np.nan
    other_pet_prev = other_pet['is_positive'].mean() * 100 if not other_pet.empty else np.nan
    
    record = {
        'Display_Element': display_element,
        'Best_Hit_ARO': aro,
        'Category': category,
        'Cut_Off': cutoff,
        'This_Study': shanghai_pet_prev,
        'Current_Pet': other_pet_prev,
        'Current_Colony': prevalence.get('Dog Colony', np.nan),
        'Current_Shelter': prevalence.get('Dog Shelter', np.nan),
        'Current_Free': prevalence.get('Dog Free_roaming', np.nan)
    }
    prevalence_records.append(record)

# Create dataframe
df = pd.DataFrame(prevalence_records)

# Extract size in kbp from df_Extra_chromo and merge with df
df_Extra_chromo['Size_kbp'] = df_Extra_chromo['Working_header'].str.extract(r'Size_(\d+)').astype(float) / 1000
size_data = df_Extra_chromo[['Element', 'Size_kbp']]
df = pd.merge(df, size_data, left_on='Display_Element', right_on='Element', how='left')
df = df.drop(columns=['Element'])

# Load and process putative hosts
putative_hosts = pd.read_csv('intermediate-outputs/external_datasets_mappings/dogs_putative_hosts.csv', sep=',')
host_mapping = putative_hosts.groupby('Putative Plasmid')['Putative host'].first().str.extract(r's__([^;]+)').iloc[:, 0].to_dict()
df['Putative_Host'] = df['Display_Element'].str.replace('EC', 'NC').map(host_mapping).fillna('-').str.replace('_', ' ')

# Aggregate by Display_Element, combining ARGs and Cut_Off values
df = df.groupby('Display_Element').agg({
    'Best_Hit_ARO': lambda x: ', '.join(x.dropna().unique()),
    'Cut_Off': lambda x: list(x.dropna()),
    'Category': 'first',
    'This_Study': 'first',
    'Current_Pet': 'first',
    'Current_Colony': 'first',
    'Current_Shelter': 'first',
    'Current_Free': 'first',
    'Size_kbp': 'first',
    'Putative_Host': 'first'
}).reset_index()

# Set up plot
n_rows = len(df)
fig_height = max(12 / IN2CM, n_rows * 0.7 / IN2CM)
fig = plt.figure(figsize=(13 / IN2CM, fig_height))
gs = gridspec.GridSpec(2, 6, height_ratios=[6, 0.5], width_ratios=[0.5, 1.1, 1.8, 0.1, 0.7, 0.9], hspace=0.3, wspace=0.001)

# colormap
prevalence_cmap = cm.YlOrBr
prevalence_cmap.set_bad(color='#DCDCDC') #lightgrey color

# Plot Category and Size
ax0 = fig.add_subplot(gs[0, 0])
ax0.set_xlim(10, 500)
ax0.set_xscale('log')
ax0.set_ylim(0, n_rows)
ax0.set_xlabel('Size (kbp)', fontsize=9, labelpad=5)
ax0.set_xticks([10, 50, 100])
ax0.tick_params(axis='x', labelsize=9)
ax0.set_xticklabels(['10', '50', '100'], fontsize=9)
ax0.set_yticks([])
ax0.invert_xaxis()
ax0.spines['left'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
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
for i, name in enumerate(df['Display_Element']):
    ax1.text(0, n_rows - i - 0.5, name, va='center', fontsize=9)
ax1.text(0.5, n_rows + 0.5, "Element", fontsize=9, ha='center', rotation=45)

# Plot Heatmap with square cells, with Non prevalent elements as light grey color
prevalence_columns = ['This_Study', 'Current_Pet', 'Current_Colony', 'Current_Shelter', 'Current_Free']
heatmap_data = df[prevalence_columns].copy()
heatmap_data.columns = ['This_Study', 'Pet', 'Colony', 'Shelter', 'Free']  # Updated column order
heatmap_plot_data = heatmap_data.replace(0, np.nan)
ax2 = fig.add_subplot(gs[0, 2])
sns.heatmap(heatmap_plot_data, ax=ax2, cmap=prevalence_cmap, vmin=0, vmax=50,
            cbar=False, xticklabels=True, yticklabels=False,
            linewidths=0, linecolor='none', square=True)
ax2.set_xticks(np.arange(5) + 0.5)
ax2.set_xticklabels(['This_Study', 'Pet', 'Colony', 'Shelter', 'Free'], rotation=45, ha='center', fontsize=9)  # Updated label order
ax2.xaxis.tick_top()
ax2.set_xlabel('')

# Plot Colorbar
cb_ax = fig.add_subplot(gs[1, 2])
plt.colorbar(ax2.collections[0], cax=cb_ax, orientation='horizontal')
cb_ax.tick_params(labelsize=9)

# Plot ARG labels, stacking ARGs vertically for all elements with multiple ARGs
ax3 = fig.add_subplot(gs[0, 4])
ax3.set_xlim(0, 1)
ax3.set_ylim(0, n_rows)
ax3.axis('off')
line_height = 0.5
for i, (element, arg_str, cutoffs) in enumerate(zip(df['Display_Element'], df['Best_Hit_ARO'], df['Cut_Off'])):
    args = [a.strip() for a in arg_str.split(',')]
    if len(args) > 1:  # Stack ARGs vertically if multiple
        for j, (arg, cutoff) in enumerate(zip(args, cutoffs)):
            weight = 'bold' if cutoff == 'Perfect' else 'normal'
            ax3.text(0, n_rows - i - 0.5 + (len(args)/2 - j - 1) * line_height, arg, va='center', fontsize=9, fontweight=weight)
    else:  # Single ARG
        weight = 'bold' if cutoffs[0] == 'Perfect' else 'normal'
        ax3.text(0, n_rows - i - 0.5, arg_str, va='center', fontsize=9, fontweight=weight)
ax3.text(0, n_rows + 0.5, "ARG", fontsize=9, ha='left', rotation=45)

# Plot Putative Host
ax4 = fig.add_subplot(gs[0, 5])
ax4.set_xlim(0, 1)
ax4.set_ylim(0, n_rows)
ax4.axis('off')
for i, host in enumerate(df['Putative_Host']):
    ax4.text(0, n_rows - i - 0.5, host, va='center', fontsize=9)
ax4.text(0, n_rows + 0.5, "Putative Host", fontsize=9, ha='left', rotation=45)

# Finalize and save plot
fig.tight_layout()
fig.subplots_adjust(top=0.95, bottom=0.10, hspace=0.3, wspace=0.001)
fig.savefig('prevalence_heatmap.svg', format='svg', bbox_inches='tight', dpi=300)

# Part 2: Circular Gene Plot 
faa_file = "intermediate-outputs/Prodigal/D003/D003_proteins.faa.gz"
annotations_file = "intermediate-outputs/eggNOG_annot_contigs/D003/D003.emapper.annotations"
fasta_file = "shanghai_dogs/data/ShanghaiDogsAssemblies/D003_PP1_PolcaCorr.fna.gz"

#Full contig length from fasta.gz
def get_full_contig_length(contig_prefix, fasta_gz_path):
    with gzip.open(fasta_gz_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id.startswith(contig_prefix):
                return len(record.seq)
    return None
contig_prefix = "contig_2645_polypolish"

# Get full contig length 
genomic_size = get_full_contig_length(contig_prefix, fasta_file)
if genomic_size is None:
    raise ValueError(f"Contig {contig_prefix} not found in {fasta_file}")

# Step 1: Parse gene coordinates from .faa.gz
pattern = re.compile(r"#\s+(\d+)\s+#\s+(\d+)")
genes = []

with gzip.open(faa_file, "rt") as f:
    for line in f:
        if line.startswith(f">{contig_prefix}"):
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
valid_cogs = set("DLISJLV")
cog_categories = {}

with open(annotations_file, "r") as f:
    for line in f:
        if line.startswith("#") or not line.startswith(contig_prefix):
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

# Step 4: Plot circular gene map
fig, ax = plt.subplots(figsize=(3, 3))  
fig_circle = plt.Circle((0, 0), 0.15, fill=False, color='black', linewidth=1.5)
ax.add_artist(fig_circle)

# Color map setup
unique_cogs = sorted(set(g["cog_category"] for g in genes))
color_map = {"S": [0.6, 0.6, 0.6]}  # grey for 'S'
non_s_cogs = [c for c in unique_cogs if c != "S"]
color_map.update(zip(non_s_cogs, [c for i, c in enumerate(plt.cm.Dark2.colors) if i != 7][:len(non_s_cogs)]))

# Plot each gene arc with strand direction
for gene in genes:
    angle = ((gene["start"] + gene["end"]) / 2 / genomic_size) * 2 * np.pi
    arc_len = (gene["end"] - gene["start"]) / genomic_size * 2 * np.pi
    theta = np.linspace(angle - arc_len / 2, angle + arc_len / 2, 100)
    r = 0.15
    x, y = r * np.cos(theta), r * np.sin(theta)
    cog = gene["cog_category"]
    color = color_map.get(cog)

    if gene["strand"] == -1:
        x, y = x[::-1], y[::-1]

    ax.plot(x, y, linewidth=4, color=color)

    strand_symbol = "+" if gene["strand"] == 1 else "-"
    label = f"{cog}{strand_symbol}" if cog != "NA" else f"NA{strand_symbol}"
    ax.text(0.19 * np.cos(angle), 0.19 * np.sin(angle), label, ha='center', va='center', fontsize=9, color='black')

ax.set_aspect('equal')
ax.set_xlim(-0.3, 0.3)
ax.set_ylim(-0.3, 0.3)
ax.axis("off")

fig.tight_layout()
fig.savefig('circular_gene_plot.svg')
