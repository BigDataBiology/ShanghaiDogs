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

# Clean contig names to match
df_arg['Contig_clean'] = df_arg['Contig'].str.extract(r'(contig_\d+)')

# Merge ARGs with extrachromosomal elements
merged_df = pd.merge(
    df_non_chromo,
    df_arg[['Contig_clean', 'sample_id', 'Best_Hit_ARO', 'Cut_Off', 'Best_Identities']],
    left_on=['Contig', 'Sample'],
    right_on=['Contig_clean', 'sample_id'],
    how='left'
)

merged_df = merged_df.drop(columns=['Contig_clean'])

# Filter rows with ARO information
filtered_df = merged_df[merged_df['Best_Hit_ARO'].notnull()]

# Build prevalence dataframe
# Create one row per (Element, ARO) pair
element_info_rows = filtered_df[['Element', 'Best_Hit_ARO']].drop_duplicates()
element_info_list = list(element_info_rows.itertuples(index=False, name=None))

abbreviations = {
    'Dog Pet': 'Current_Pet',
    'Dog Colony': 'Current_Colony',
    'Dog Shelter': 'Current_Shelter',
    'Dog Free_roaming': 'Current_Free'
}

prevalence_records = []

for element_id, aro in element_info_list:
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
        'Current_Pet': prevalence.get('Dog Pet', np.nan),
        'Current_Colony': prevalence.get('Dog Colony', np.nan),
        'Current_Shelter': prevalence.get('Dog Shelter', np.nan),
        'Current_Free': prevalence.get('Dog Free_roaming', np.nan)
    }
    prevalence_records.append(record)

df = pd.DataFrame(prevalence_records)

# Plotting 
fig = plt.figure(figsize=(9 / IN2CM, 7 / IN2CM))  # Adjusted to fit ~9 cm width
gs = gridspec.GridSpec(2, 2, height_ratios=[6, 0.3], width_ratios=[1, 4])

# Colormap
ylorbr_colors = cm.YlOrBr(np.linspace(0, 1, 256))  # Smooth gradient
prevalence_cmap = LinearSegmentedColormap.from_list('YlOrBr_Custom', ylorbr_colors)

#  ARO labels column
ax1 = plt.subplot(gs[0, 0])
ax1.imshow(np.ones((len(df), 1)), cmap=LinearSegmentedColormap.from_list("white_only", ['white', 'white']),
           aspect='auto')  # aspect='auto' to minimize gaps
ax1.set_yticks(np.arange(len(df)))
ax1.set_yticklabels(df['Element'], fontsize=9)
ax1.set_xticks([])
ax1.set_title('ARG', fontsize=9, rotation=45, ha='right')
ax1.xaxis.tick_top()
for i, aro in enumerate(df['Best_Hit_ARO']):
    ax1.text(0, i, aro, ha='center', va='center', fontsize=9, color='black')
for spine in ax1.spines.values():
    spine.set_visible(False)

# Single Prevalence Heatmap using Seaborn 
prevalence_columns = ['Current_Pet', 'Current_Colony', 'Current_Free', 'Current_Shelter']
full_titles = ['Pet', 'Colony', 'Free-roaming', 'Shelter']
heatmap_data = df[prevalence_columns].copy()
heatmap_data.columns = full_titles

ax2 = plt.subplot(gs[0, 1])
sns.heatmap(heatmap_data, ax=ax2, cmap=prevalence_cmap, vmin=0, vmax=50, cbar=False,
            annot=False, linewidths=0, linecolor='none')
ax2.set_yticks([])
ax2.set_yticklabels([])
ax2.set_xticklabels(full_titles, rotation=45, ha='left', fontsize=9)
ax2.xaxis.tick_top()
ax2.set_xlabel('')
for spine in ax2.spines.values():
    spine.set_visible(False)

# colorbar
cb_ax = plt.subplot(gs[1, 1])
cb = plt.colorbar(ax2.collections[0], cax=cb_ax, orientation='horizontal')
cb.set_label('Prevalence (%)', fontsize=9)
cb.ax.tick_params(labelsize=9)

# Layout adjustments
fig.tight_layout()
fig.subplots_adjust(top=0.92, bottom=0.12, hspace=0.3, wspace=0.05)
fig.show()

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
                contig_id = line.split()[0][1:]  # remove ">"
                genes.append({
                    "start": start,
                    "end": end,
                    "contig": contig_id
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
            cog_categories[contig_name] = filtered if filtered else "NA"

# Step 3: Map COGs to genes
for gene in genes:
    gene["cog_category"] = cog_categories.get(gene["contig"], "NA")

# Step 4: Plot
genomic_size = max(g["end"] for g in genes)
fig, ax = plt.subplots(figsize=(2, 2))  # Size remains reduced by 1.5 times

fig_circle = plt.Circle((0, 0), 0.15, fill=False, color='black')
ax.add_artist(fig_circle)

# Color map: Dark2 greyish-green for NA, Dark2 for others
unique_cogs = sorted(set(g["cog_category"] for g in genes))
color_map = {"NA": plt.cm.Dark2.colors[7]}
non_na_cogs = [c for c in unique_cogs if c != "NA"]
color_map.update(zip(non_na_cogs, [c for i, c in enumerate(plt.cm.Dark2.colors) if i != 7][:len(non_na_cogs)]))

# Arcs and labels
for gene in genes:
    angle = ((gene["start"] + gene["end"]) / 2 / genomic_size) * 2 * np.pi
    arc_len = (gene["end"] - gene["start"]) / genomic_size * 2 * np.pi
    theta = np.linspace(angle - arc_len / 2, angle + arc_len / 2, 100)
    x, y = 0.15 * np.cos(theta), 0.15 * np.sin(theta)
    cog = gene["cog_category"]
    color = color_map.get(cog)
    ax.plot(x, y, linewidth=4, color=color)
    label = "NA" if cog == "NA" else cog
    ax.text(0.19 * np.cos(angle), 0.19 * np.sin(angle), label, ha='center', va='center', fontsize=7)

# Formatting
ax.set_aspect('equal')
ax.set_xlim(-0.25, 0.25)  
ax.set_ylim(-0.25, 0.25) 
ax.axis("off")

# Legend
handles = [plt.Line2D([0], [0], color=color_map[c], lw=4, label=c) for c in unique_cogs]
fig.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=9)

fig.tight_layout()
fig.show()
