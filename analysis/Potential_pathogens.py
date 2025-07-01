import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
import os
import argnorm.lib
import argnorm.drug_categorization

def antibiotic(aro):
    if pd.isna(aro):
        return 'Unknown'
    aro = f'ARO:{int(aro)}'
    drugs = argnorm.drug_categorization.confers_resistance_to(aro)
    drug_classes = argnorm.drug_categorization.drugs_to_drug_classes(drugs)
    names = [aro_ontology[d].name for d in drug_classes]
    names = set(names)
    names = [n.replace(' antibiotic', '') for n in names]
    names.sort()
    return ';'.join(names)

def extract_species(classification):
    if pd.isna(classification):
        return "Unknown"
    genus, species = None, None
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

mags = pd.read_csv("data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv")
args_mags = pd.read_csv("intermediate-outputs/06_ARG/MAGs-ARGs_ALL_filt.txt")

merged_df = pd.merge(mags, args_mags, on="Bin ID")
merged_df['Species'] = merged_df['Classification'].map(extract_species)
merged_df['Species'].unique()

species = merged_df['Species'].to_list()
potential_pathogens = [sp for sp in species if sp.startswith('Helicobacter')] + \
                     [sp for sp in species if sp.startswith('Enterococcus')] + \
                     [sp for sp in species if sp.startswith('Staphylococcus')] + \
                     ['Escherichia coli', 'Proteus mirabilis', 'Clostridioides difficile', 'Sarcina ventriculi', 'Klebsiella pneumoniae']
merged_df['Potential_Pathogen'] = merged_df['Species'].apply(set(potential_pathogens).__contains__)

merged_df = merged_df[~merged_df['ARO'].isna()]
merged_df = merged_df.query('Potential_Pathogen')
aro_ontology = argnorm.lib.get_aro_ontology()
resf = argnorm.lib.get_aro_mapping_table('resfinder')
in_resfinder = set(resf['ARO'].str.split(':').str[1].map(float))
in_resfinder.discard(float('nan'))
merged_df = merged_df[merged_df['ARO'].map(set(in_resfinder).__contains__)]
merged_df = merged_df.query('Best_Identities >= 95')

merged_df['antibiotic'] = merged_df['ARO'].map(antibiotic)
merged_df['Gene_name'] = merged_df['Best_Hit_ARO'] + ' (' + merged_df['antibiotic'] + ')'

species_counts = merged_df['Species'].value_counts().to_dict()
merged_df['Species_Prevalence'] = merged_df['Species'].map(species_counts)

# ARG heatmap by antibiotic class
ab_heatmap = pd.crosstab(merged_df['Bin ID'], merged_df['antibiotic'])
bin2species = merged_df.set_index('Bin ID')['Species'].to_dict()
ab_heatmap.index = ab_heatmap.index.map(bin2species)
ab_heatmap = ab_heatmap.groupby(ab_heatmap.index).sum()

# prevalence column 
ab_heatmap['Prevalence'] = ab_heatmap.index.map(species_counts)

#MDR resistance to 3+ classes
ab_heatmap['Num_Drug_Classes'] = (ab_heatmap.iloc[:, :-1] > 0).sum(axis=1)
ab_heatmap['MDR'] = ab_heatmap['Num_Drug_Classes'] >= 3
ab_heatmap.to_csv('ARG_species_summary.csv')

# Plotheatmap
mheat = ab_heatmap.iloc[:, :-3]
mheat[mheat == 0] = np.nan
fig, ax = plt.subplots(figsize=(17 / 2.54, 10 / 2.54))
cm.OrRd.set_bad(color='white')
sns.heatmap(mheat, cmap='OrRd', cbar_kws={'label': 'ARG Count'}, linewidths=0.5, ax=ax)
ax.set_xlabel('')
ax.set_ylabel('')
fig.tight_layout()
fig.savefig('Updated_hetamap.png', dpi=300)
