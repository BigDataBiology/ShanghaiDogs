import numpy as np
from matplotlib import cm
import matplotlib.patches as mpatches
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argnorm.lib
import argnorm.drug_categorization
from collections import Counter

mags = pd.read_csv("data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv")
args_mags = pd.read_csv("intermediate-outputs/06_ARG/MAGs-ARGs_ALL_filt.txt")
merged_df = pd.merge(mags, args_mags, on="Bin ID")

def extract_species(classification):
    if pd.isna(classification):
        return "Unknown"
    genus, species = "", ""
    for part in classification.split(';'):
        if part.startswith('g__'):
            genus = part[3:]
        if part.startswith('s__'):
            species = part[3:]
    return species if species else f"{genus} novel_sp" if genus else "Unknown"

merged_df['Species'] = merged_df['Classification'].map(extract_species)

potential_pathogens = [
    'Escherichia coli', 'Proteus mirabilis', 'Clostridioides difficile',
    'Sarcina ventriculi', 'Klebsiella pneumoniae'
] + [s for s in merged_df['Species'].unique()
     if s.startswith(('Helicobacter', 'Enterococcus', 'Staphylococcus'))]

mag_counts_df = (
    merged_df[merged_df['Species'].isin(potential_pathogens)][['Species', 'Bin ID']]
    .drop_duplicates()
    .groupby('Species')
    .count()
    .rename(columns={'Bin ID': 'MAG Count'})
)
mag_counts = mag_counts_df['MAG Count'].to_dict()

merged_df = merged_df[
    merged_df['Species'].isin(potential_pathogens) &
    (~merged_df['ARO'].isna())
]

aro_ontology = argnorm.lib.get_aro_ontology()
resf = argnorm.lib.get_aro_mapping_table('resfinder')
in_resfinder = set(resf['ARO'].str.split(':').str[1].map(float))
merged_df = merged_df[
    merged_df['ARO'].map(set(in_resfinder).__contains__) &
    (merged_df['Best_Identities'] >= 95)
]

def antibiotic(aro):
    if pd.isna(aro):
        return 'Unknown'
    aro_id = f'ARO:{int(aro)}'
    drugs = argnorm.drug_categorization.confers_resistance_to(aro_id)
    drug_classes = argnorm.drug_categorization.drugs_to_drug_classes(drugs)
    names = [aro_ontology[d].name.replace(' antibiotic', '') for d in drug_classes]
    return ';'.join(sorted(set(names)))

merged_df['antibiotic'] = merged_df['ARO'].map(antibiotic)

def abbreviate_species(name):
    parts = name.split()
    if len(parts) >= 2:
        return f"$\\it{{{parts[0][0]}. {' '.join(parts[1:])}}}$"
    return name

merged_df['SpeciesWithCounts'] = merged_df['Species'].map(
    lambda s: f"{abbreviate_species(s)} – {mag_counts.get(s, 0)}"
)

heatmap = pd.crosstab(merged_df['Bin ID'], merged_df['Best_Hit_ARO'])
bin2species = merged_df.set_index('Bin ID')['SpeciesWithCounts'].to_dict()
heatmap.index = heatmap.index.map(bin2species)
heatmap.sort_index(inplace=True)

arg_to_class = merged_df.set_index('Best_Hit_ARO')['antibiotic'].to_dict()
arg_primary_class = {
    arg: arg_to_class.get(arg, 'Unknown').split(';')[0]
    if arg_to_class.get(arg) else 'Unknown'
    for arg in heatmap.columns
}
class_counts = Counter(arg_primary_class.values())
sorted_classes = [cls for cls, _ in class_counts.most_common()]
sorted_args = sorted(
    heatmap.columns,
    key=lambda x: (sorted_classes.index(arg_primary_class[x]), x)
)
heatmap = heatmap[sorted_args]

mheat = heatmap.groupby(heatmap.index).mean()
mheat[mheat == 0] = np.nan

species_order = sorted(
    mheat.index,
    key=lambda s: int(s.split(' – ')[-1]),
    reverse=True
)
mheat = mheat.loc[species_order]

def format_species_label(species_label):
    try:
        name_part, count = species_label.split(' – ')
    except ValueError:
        return species_label
    return f"{name_part} (n={count})"

mheat.index = mheat.index.map(format_species_label)

def italicize_gene_name(gene):
    if "AAC(6')-Ie-APH(2'')-Ia bifunctional protein" in gene:
        gene = gene.replace(" bifunctional protein", "")
    return r"$\it{" + gene.replace("_", r"\_") + "}$"

mheat.columns = [italicize_gene_name(col) for col in mheat.columns]

all_classes = sorted(set(arg_primary_class.values()), key=lambda x: sorted_classes.index(x))
palette = sns.color_palette('Dark2', len(all_classes))
class_to_color = dict(zip(all_classes, palette))
arg_colors = [class_to_color[arg_primary_class[arg]] for arg in sorted_args]

IN2CM = 2.54
plt.rcParams.update({'font.size': 8})
fig, ax = plt.subplots(figsize=(17 / IN2CM, 10 / IN2CM))
cmap = cm.OrRd
cmap.set_bad(color='white')
heatmap_img = ax.imshow(mheat, aspect='equal', cmap=cmap, interpolation='nearest')

ax.spines[:].set_visible(False)
ax.set_xticks(np.arange(len(mheat.columns) + 1) - 0.5, minor=True)
ax.set_yticks(np.arange(len(mheat.index) + 1) - 0.5, minor=True)
ax.grid(which="minor", color="w", linestyle='-', linewidth=0.5)
ax.tick_params(which="minor", bottom=False, left=False)
ax.set_xticks(range(len(mheat.columns)))
ax.set_xticklabels(mheat.columns, rotation=90, fontsize=8)
ax.set_yticks(range(len(mheat.index)))
ax.set_yticklabels(mheat.index, fontsize=8)

for idx, color in enumerate(arg_colors):
    ax.add_patch(plt.Rectangle((idx - 0.5, len(mheat.index) - 0.4), 1, 0.3,
                               transform=ax.transData, clip_on=False,
                               facecolor=color, linewidth=0))

handles = [mpatches.Patch(color=class_to_color[c], label=c) for c in class_to_color]
legend = ax.legend(
    handles=handles,
    loc='upper center',
    bbox_to_anchor=(0.8, -0.26),
    ncol=2,
    frameon=False,
    fontsize=8,
    title='Drug Class',
    title_fontsize=8
)

cbar = fig.colorbar(heatmap_img, ax=ax, shrink=0.3, aspect=18, pad=0.02)
cbar.ax.tick_params(labelsize=8)
sns.despine(ax=ax, trim=True)
ax.tick_params(axis='x', which='major', length=0)
fig.tight_layout()
fig.savefig('potential_pathogens_updated.svg', dpi=300, bbox_inches='tight')
