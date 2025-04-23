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

#constant
threshold = 0.8

# Load metadata
metadata = pd.read_csv('data/ShanghaiDogsMetadata/REP_canid_metadata.csv', sep=',', encoding='unicode_escape', index_col=0)
metadata['Sex'] = metadata['Sex'].str.lower()

# Load coverage and NCE info
MAGs_NCE_covered_frac = pd.read_csv('intermediate-outputs/external_datasets_mappings/SHD_covered_fraction.tsv.gz',
                                    sep='\t', index_col=0)
NCE_info = pd.read_csv('data/ShanghaiDogsTables/SHD_NC_props.tsv.gz',
                       sep='\t', index_col=0)

# --- Define elements and best-hit AROs ---
element_info = {
    'SHD1_NC.006': 'OXA-85',
    'SHD1_NC.021': 'OXA-85',
    'SHD1_NC.026': 'TEM-1',
    'SHD1_NC.053': 'ErmQ',
    'SHD1_NC.094': 'lnuC',
    'SHD1_NC.113': 'lnuC',
    'SHD1_NC.143': 'lnuC',
    'SHD1_NC.147': 'lnuC',
    'SHD1_NC.174': 'tet(Q)',
    'SHD1_NC.183': 'ErmQ / tetB(P) / tetA(P)'
}
element_ids = list(element_info.keys())

abbreviations = {
    'Dog Pet': 'Current_Pet',
    'Dog Colony': 'Current_Colony',
    'Dog Shelter': 'Current_Shelter',
    'Dog Free_roaming': 'Current_Free'
}

# --- Compute prevalence ---
prevalence_records = []

for element_id in element_ids:
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
        'Best_Hit_ARO': element_info[element_id],
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
ylorrd_colors = cm.get_cmap('YlOrRd', 8)(range(8))
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
    "Prevalence threshold: ≥ 80%"
)
ax_right.text(0.2, 1, legend_text, fontsize=6, va='top', ha='left', linespacing=1.5)

fig.tight_layout()
fig.subplots_adjust(top=0.88, wspace=0.4)
fig.show()

# --- Part 2: Circular Gene Plot ---
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
