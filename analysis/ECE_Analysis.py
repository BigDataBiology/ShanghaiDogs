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

# matplotlib uses inches, so we need to convert cm to inches
IN2CM = 2.54

# --- Part 1: Prevalence Bar Plot ---
# Load metadata
metadata = pd.read_csv('data/ShanghaiDogsMetadata/REP_canid_metadata.csv', sep=',', encoding='unicode_escape', index_col=0)
metadata['Sex'] = metadata['Sex'].str.lower()

# Load coverage and NCE info
MAGs_NCE_covered_frac = pd.read_csv('intermediate-outputs/external_datasets_mappings/SHD_covered_fraction.tsv.gz',
                                    sep='\t', index_col=0)
NCE_info = pd.read_csv('data/ShanghaiDogsTables/SHD_NC_props.tsv.gz',
                       sep='\t', index_col=0)

# Define your element and threshold
element_id = 'SHD1_NC.006'
threshold = 0.8

# Abbreviations for plotting
abbreviations = {
    'Dog Pet': 'P',
    'Dog Colony': 'C',
    'Dog Shelter': 'S',
    'Dog Free_roaming': 'F'
}

assert (element_id in NCE_info.index and element_id in MAGs_NCE_covered_frac.index), \
    f"Element {element_id} not found in NCE info or coverage data."

# Get coverage and apply threshold
element_coverage = MAGs_NCE_covered_frac.loc[element_id]
element_presence = element_coverage[element_coverage >= threshold]
positive_samples = element_presence.index.str.split('_').str[0]

# Get metadata for all and positive samples
metadata_env = metadata[['env_classification']].copy()
metadata_env = metadata_env[metadata_env['env_classification'].isin(abbreviations.keys())]
metadata_env['is_positive'] = metadata_env.index.isin(positive_samples)

# Compute prevalence
prevalence = metadata_env.groupby('env_classification')['is_positive'].mean().reindex(abbreviations.keys(), fill_value=0) * 100
counts = metadata_env.groupby('env_classification')['is_positive'].sum().reindex(abbreviations.keys(), fill_value=0)

# Print prevalence
print("Prevalence of SHD1_NC.006 (≥0.8) per environment:")
for env, pct in prevalence.items():
    print(f"  {env:16} ({abbreviations[env]}): {pct:.1f}%")

# Plot
fig_width = 12 / IN2CM
fig_height = 6 / IN2CM

fig, ax = plt.subplots(figsize=(fig_width, fig_height))  # Set figsize in inches
x_labels = [abbreviations[env] for env in prevalence.index]
palette = sns.color_palette('Dark2', n_colors=len(prevalence))
sns.barplot(x=x_labels, y=prevalence.values, palette=palette, ax=ax)

# % labels
for i, pct in enumerate(prevalence.values):
    ax.text(i, pct + 1, f'{pct:.1f}%', ha='center', va='bottom', fontsize=10)

ax.set_label('Prevalence (%)')
ax.set_label('Environment')
ax.set_ylim(0, max(prevalence.values) * 1.2 + 5)

# Colored legend
legend_elements = [
    Patch(facecolor=palette[i], edgecolor='black', label=f"{abbreviations[env]} = {env}")
    for i, env in enumerate(prevalence.index)
]
ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

fig.tight_layout()
fig.show()


# --- Part 2: Heatmap Visualization ---

# Data for heatmap
data = {
    'Element': ['SHD1_NC.006', 'SHD1_NC.021', 'SHD1_NC.026', 'SHD1_NC.053'],
    'Best_Hit_ARO': ['OXA-85', 'OXA-85', 'TEM-1', 'ErmQ'],
    'Current_Pet': [18.8, np.nan, np.nan, np.nan],
    'Current_Colony': [14.7, np.nan, np.nan, np.nan],
    'Current_Shelter': [0.0, np.nan, np.nan, np.nan],
    'Current_Free': [2.8, np.nan, np.nan, np.nan]
}
df = pd.DataFrame(data)

# Set figure size
fig = plt.figure(figsize=(12./IN2CM, 6/IN2CM))  # ~12cm x 6cm
gs = gridspec.GridSpec(1, 2, width_ratios=[3.5, 1])

# Subgrid for heatmap
gs_heatmap = gridspec.GridSpecFromSubplotSpec(1, 5, subplot_spec=gs[0], width_ratios=[1, 1, 1, 1, 1])

# Colors from Dark2 colormap
dark2_colors = cm.get_cmap('Dark2', 8)(range(8))
aro_color = dark2_colors[1]
prevalence_cmap = LinearSegmentedColormap.from_list('Dark2_Prev', dark2_colors)

# ARO column
ax1 = plt.subplot(gs_heatmap[0])
aro_mapping = {aro: 1.0 if aro == 'OXA-85' else 0.0 for aro in df['Best_Hit_ARO']}
aro_numeric = np.array([aro_mapping[aro] for aro in df['Best_Hit_ARO']]).reshape(-1, 1)
im1 = ax1.imshow(aro_numeric, cmap=LinearSegmentedColormap.from_list('ARO_cmap', [aro_color, dark2_colors[3]]), aspect=0.4)
ax1.set_yticks(np.arange(len(df)))
ax1.set_yticklabels(df['Element'], fontsize=5)
ax1.set_xticks([])
ax1.set_title('ARO', fontsize=6)
ax1.xaxis.tick_top()
for i, aro in enumerate(df['Best_Hit_ARO']):
    ax1.text(0, i, aro, ha='center', va='center', fontsize=5, weight='normal')

# Prevalence columns for current study only
prevalence_columns = ['Current_Pet', 'Current_Colony', 'Current_Free', 'Current_Shelter']
simple_titles = ['P', 'C', 'F', 'S']

# Create new axes for each prevalence column
for i, (col, title) in enumerate(zip(prevalence_columns, simple_titles)):
    ax = plt.subplot(gs_heatmap[i+1])  # Create new subplot for each prevalence column

    # Mask 0.0 values
    prevalence_data = df[col].replace(0.0, np.nan).values.reshape(-1, 1)  # Replace 0.0 with NaN
    masked = np.ma.masked_invalid(prevalence_data)  # Mask NaN values

    im = ax.imshow(masked, cmap=prevalence_cmap, aspect=0.4, vmin=0, vmax=50)
    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_title(title, fontsize=6)
    ax.xaxis.tick_top()

# Right-side legend panel
ax_right = plt.subplot(gs[1])
ax_right.axis('off')
legend_text = (
    "ARO: Best hit antimicrobial resistance gene\n"
    "P: Dog Pet\n"
    "C: Dog Colony\n"
    "F: Dog Free-roaming\n"
    "S: Dog Shelter"
)

ax_right.text(0, 1, legend_text, fontsize=6, va='top', ha='left', linespacing=1.5)

# Layout adjustments
fig.tight_layout()
fig.subplots_adjust(top=0.88, wspace=0.3)
fig.savefig("heatmap_single_legend.svg", format='svg')
fig.show()

# --- Part 3: Circular Gene Plot ---
faa_file = "intermediate-outputs/Prodigal/D003/D003_proteins.faa.gz"

# Initialize gene list
genes = []

# Pattern to extract coordinates from headers like:
# >contig_2645_polypolish_1 # 202 # 690 # 1 ...
pattern = re.compile(r"#\s+(\d+)\s+#\s+(\d+)")

# Read .faa.gz file and extract gene start/end positions for contig_2645
with gzip.open(faa_file, "rt") as f:
    for line in f:
        if line.startswith(">contig_2645"):
            match = pattern.search(line)
            if match:
                start, end = map(int, match.groups())
                genes.append({"start": start, "end": end})

# Check if any genes were found
if not genes:
    raise ValueError("No gene entries found for contig_2645.")

# Calculate total genome size (for angle scaling)
genomic_size = max(g["end"] for g in genes)

# Create circular plot
fig, ax = plt.subplots(figsize=(4, 4))

# Draw main circle
circle = plt.Circle((0, 0), 1, fill=False, color='black')
ax.add_artist(circle)

# Color map for genes
colors = plt.cm.Dark2(np.linspace(0, 1, len(genes)))

# Draw arcs for each gene
for i, (gene, color) in enumerate(zip(genes, colors)):
    mean_pos = (gene["start"] + gene["end"]) / 2
    angle = (mean_pos / genomic_size) * 2 * np.pi
    length = (gene["end"] - gene["start"]) / genomic_size * 2 * np.pi
    theta = np.linspace(angle - length/2, angle + length/2, 100)
    x = np.cos(theta)
    y = np.sin(theta)
    ax.plot(x, y, linewidth=5, color=color, label=f'Gene {i+1}')

    # Add gene label outside the circle
    label_x = 1.1 * np.cos(angle)
    label_y = 1.1 * np.sin(angle)
    ax.text(label_x, label_y, f'G{i+1}', ha='center', va='center', fontsize=8)

# Style adjustments
ax.set_aspect('equal')
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.axis('off')

# Show legend
fig.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=6)
fig.tight_layout()

# Show plot (or replace with savefig if needed)
fig.show()
