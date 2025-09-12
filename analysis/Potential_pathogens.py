import numpy as np
from matplotlib import cm
import matplotlib.patches as mpatches
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argnorm.lib
import argnorm.drug_categorization
from collections import Counter
plt.rcParams['svg.fonttype'] = 'none'


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


mags = pd.read_csv("data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv")
args_mags = pd.read_csv("intermediate-outputs/06_ARG/MAGs-ARGs_ALL_filt.txt")
mags['Species'] = mags['Classification'].map(extract_species)

merged_df = pd.merge(mags, args_mags, on="Bin ID", how="left")

potential_pathogens = [
    'Klebsiella pneumoniae',
    'Proteus mirabilis',
    'Escherichia coli',
    'Clostridioides difficile',
    'Citrobacter freundii',
    'Citrobacter portucalensis',
    'Streptococcus lutetiensis',
] + [s for s in merged_df['Species'].unique()
     if s.startswith(('Helicobacter', 'Enterococcus', 'Staphylococcus', 'Campylobacter'))]

potential_pathogens = set(potential_pathogens)

species_counts = (
     mags[mags['Species'].isin(potential_pathogens)][['Species', 'Bin ID']]
    .drop_duplicates()
    .groupby('Species')
    .count()
)['Bin ID'].to_dict()

merged_df = merged_df[
    merged_df['Species'].isin(potential_pathogens) &
    (~merged_df['ARO'].isna())
]

#Resfinder mapping
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
    names = [aro_ontology[d].name.replace(' antibiotic', '').replace(' ; ', ';').replace(' ', '') for d in drug_classes]
    return ';'.join(sorted(set(names)))

def categorize_antibiotic_class(abx_str):
    if abx_str == 'Unknown' or pd.isna(abx_str):
        return 'Unknown'
    if ';' in abx_str:
        return 'Multi-class'
    return abx_str

merged_df['antibiotic'] = merged_df['ARO'].map(antibiotic)
merged_df['antibiotic_class_cat'] = merged_df['antibiotic'].map(categorize_antibiotic_class)


# generate ARG-by-MAG heatmap table
heatmap = pd.crosstab(merged_df['Bin ID'], merged_df['Best_Hit_ARO'])

bin2species = merged_df.set_index('Bin ID')['Species'].to_dict()
heatmap.index = heatmap.index.map(bin2species)
heatmap.sort_index(inplace=True)

arg_to_class = merged_df.set_index('Best_Hit_ARO')['antibiotic_class_cat'].to_dict()
arg_primary_class = {arg: arg_to_class.get(arg, 'Unknown') for arg in heatmap.columns}
class_counts = Counter(arg_primary_class.values())
sorted_classes = [cls for cls, _ in class_counts.most_common()]
sorted_args = sorted(
    heatmap.columns,
    key=lambda x: (sorted_classes.index(arg_primary_class[x]), x)
)
heatmap = heatmap[sorted_args]
heatmap[heatmap == 0] = np.nan

heatmap.columns = heatmap.columns.str.replace(" bifunctional protein", "")
heatmap.index = heatmap.index.str.replace('_', ' ')

#Generate color palette
all_classes = sorted(set(arg_primary_class.values()), key=lambda x: sorted_classes.index(x))
dark2_colors = sns.color_palette("Dark2", 8)
tab10_colors = sns.color_palette("tab10", len(all_classes) - len(dark2_colors))
combined_palette = dark2_colors + tab10_colors

#separate color for multi class drugs
if 'Multi-class' in all_classes:
    multi_class_color = (0.5, 0.5, 0.5)
    class_to_color_temp = dict(zip(all_classes, combined_palette[:len(all_classes)]))
    class_to_color_temp['Multi-class'] = multi_class_color
    class_to_color = class_to_color_temp
else:
    class_to_color = dict(zip(all_classes, combined_palette[:len(all_classes)]))

arg_colors = [class_to_color[arg_primary_class[arg]] for arg in sorted_args]


IN2CM = 2.54
plt.rcParams.update({'font.size': 7})
fig, ax = plt.subplots(figsize=(17 / IN2CM, 26 / IN2CM))
# black and white colormap
cmap = cm.get_cmap('Greys', 2)
heatmap_img = ax.imshow(heatmap, aspect='equal', cmap=cmap, interpolation='nearest', vmin=0, vmax=1)

# Format axes and grid lines
ax.spines[:].set_visible(True)
ax.spines['top'].set_visible(False)
ax.set_xticks(np.arange(len(heatmap.columns) + 1) - 0.5, minor=True)
ax.set_yticks(np.arange(len(heatmap.index) + 1) - 0.5, minor=True)
ax.grid(which="minor", color="lightgrey", linestyle='-', linewidth=0.5)
ax.tick_params(which="minor", bottom=False, left=False)
ax.set_xticks(range(len(heatmap.columns)))
ax.set_xticklabels(heatmap.columns, rotation=90, fontsize=8)
ax.tick_params(axis='x', pad=6)
ax.set_yticks(range(len(heatmap.index)))
ax.set_yticklabels(heatmap.index, fontsize=8)

for idx, color in enumerate(arg_colors):
    ax.add_patch(plt.Rectangle((idx - 0.5, -1.0), 1, 0.3,
                               transform=ax.transData, clip_on=False,
                               facecolor=color, linewidth=0))

# Legend for antibiotic classes
handles = [mpatches.Patch(color=class_to_color[c], label=c) for c in class_to_color]
legend = ax.legend(
    handles=handles,
    loc='upper center',
    bbox_to_anchor=(0.6, -0.52),
    ncol=2,
    frameon=False,
    fontsize=8,
    title_fontsize=8
)


ax.set_xlim(-0.5, len(heatmap.columns) - 0.5)
ax.set_ylim(-0.5, len(heatmap.index) - 0.5)
ax.tick_params(axis='x', which='major', length=0)

plt.subplots_adjust(bottom=0.25)
fig.tight_layout()
fig.savefig('clinically_significant.svg', dpi=300, bbox_inches='tight')
