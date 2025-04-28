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

# --- Plotting ---
fig = plt.figure(figsize=(13 / IN2CM, 6 / IN2CM))
gs = gridspec.GridSpec(1, 3, width_ratios=[4, 0.2, 1])
gs_heatmap = gridspec.GridSpecFromSubplotSpec(1, 5, subplot_spec=gs[0], width_ratios=[1, 1, 1, 1, 1])

# Colormaps
ylorrd_colors = plt.cm.Dark2(range(8))
prevalence_cmap = LinearSegmentedColormap.from_list('YlOrRd_Custom', ['white', ylorrd_colors[4]])

# ARO column
ax1 = plt.subplot(gs_heatmap[0])
ax1.imshow(np.ones((len(df), 1)), cmap=LinearSegmentedColormap.from_list("white_only", ['white', 'white']), aspect=0.4)
ax1.set_yticks(np.arange(len(df)))
ax1.set_yticklabels(df['Element'], fontsize=5)
ax1.set_xticks([])
ax1.set_title('ARO', fontsize=6)
ax1.xaxis.tick_top()
for i, aro in enumerate(df['Best_Hit_ARO']):
    ax1.text(0, i, aro, ha='center', va='center', fontsize=5, color='black')

# Prevalence columns
prevalence_columns = ['Current_Pet', 'Current_Colony', 'Current_Free', 'Current_Shelter']
simple_titles = ['P', 'C', 'F', 'S']

for i, (col, title) in enumerate(zip(prevalence_columns, simple_titles)):
    ax = plt.subplot(gs_heatmap[i + 1])
    data = df[col].replace(0.0, np.nan).values.reshape(-1, 1)
    masked = np.ma.masked_invalid(data)
    im = ax.imshow(masked, cmap=prevalence_cmap, aspect=0.4, vmin=0, vmax=50)
    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels([])  # Already shown
    ax.set_xticks([])
    ax.set_title(title, fontsize=6)
    ax.xaxis.tick_top()

# Colorbar
ax_cb = plt.subplot(gs[1])
cb = plt.colorbar(im, cax=ax_cb)
cb.set_label('Prevalence (%)', fontsize=6)
cb.ax.tick_params(labelsize=6)

# Legend panel
ax_right = plt.subplot(gs[2])
ax_right.axis('off')
legend_text = (
    "ARO: Best hit antimicrobial resistance gene\n"
    "P: Dog Pet\n"
    "C: Dog Colony\n"
    "F: Dog Free-roaming\n"
    "S: Dog Shelter\n"
)
ax_right.text(0.2, 1, legend_text, fontsize=6, va='top', ha='left', linespacing=1.5)

fig.tight_layout()
fig.subplots_adjust(top=0.88, wspace=0.4)
plt.show()

# --- Part 2: Circular Gene Plot ---
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
            cog_categories[contig_name] = filtered if filtered else "-"

# Step 3: Map COGs to genes
for gene in genes:
    gene["cog_category"] = cog_categories.get(gene["contig"], "-")

# Step 4: Plot
genomic_size = max(g["end"] for g in genes)
fig, ax = plt.subplots(figsize=(3, 3))

# Draw main circle
fig_circle = plt.Circle((0, 0), 1, fill=False, color='black')
ax.add_artist(fig_circle)

# Color map
unique_cogs = sorted(set(g["cog_category"] for g in genes if g["cog_category"] != "-"))
color_map = dict(zip(unique_cogs, plt.cm.tab10.colors[:len(unique_cogs)]))

# Arcs and labels
for gene in genes:
    angle = ((gene["start"] + gene["end"]) / 2 / genomic_size) * 2 * np.pi
    arc_len = (gene["end"] - gene["start"]) / genomic_size * 2 * np.pi
    theta = np.linspace(angle - arc_len / 2, angle + arc_len / 2, 100)
    x, y = np.cos(theta), np.sin(theta)
    cog = gene["cog_category"]
    color = color_map.get(cog, "gray")
    ax.plot(x, y, linewidth=4, color=color)
    ax.text(1.1 * np.cos(angle), 1.1 * np.sin(angle), cog, ha='center', va='center', fontsize=8)

# Final formatting
ax.set_aspect('equal')
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.axis("off")

# Legend
handles = [plt.Line2D([0], [0], color=color_map[c], lw=4, label=c) for c in unique_cogs]
fig.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)

fig.tight_layout()
fig.show()
