import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import shapiro,wilcoxon
from statsmodels.stats.multitest import multipletests

plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

### Import data
MIMAG_report = pd.read_csv('data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv', \
                    delimiter=',', header=0)
GTDB_qual = pd.read_csv('external-data/data/NCBI_genomes_ref/NCBI_genomes_qual_MIMAG_report.csv',sep=',')

# Filter Representative genomes to only include high-quality
SHD_qual_rep = MIMAG_report.query('Representative=="Yes" and Quality=="high-quality"')
merged = pd.merge(SHD_qual_rep,GTDB_qual,left_on='GTDBtk fastani Ref',right_on='Name') # include only those that are shared
merged_hq_only = merged[merged['Quality_y'].str.contains('high-quality')]

# Classify GTDBtk representative genomes
merged_hq_only['ref_new'] = merged_hq_only['GTDBtk fastani Ref']
merged_hq_only['ref_new'] = merged_hq_only['ref_new'].fillna('Novel species')

for n in merged_hq_only.index:
    ref_value = merged_hq_only.loc[n, 'ref_new']
    if 'GCF' in ref_value:
        merged_hq_only.loc[n, 'ref_new'] = 'RefSeq'
    elif 'GCA' in ref_value:
        merged_hq_only.loc[n, 'ref_new'] = 'GenBank'

merged_hq_only.to_csv("intermediate-outputs/tables/SHDvsREF-comparison-HQ.csv", index=False)

# Filter data for GenBank and RefSeq
genbank_data = merged_hq_only[merged_hq_only['ref_new'] == 'GenBank']
refseq_data = merged_hq_only[merged_hq_only['ref_new'] == 'RefSeq']

### BOXPLOTS OF CONTIGS and tRNAs: MAGs vs REF (GCA & GCF)
## Contiguity
df_ctgs = merged_hq_only[['Bin ID','Nr contigs','Number','ref_new']]
df_ctgs.columns = ['Bin ID','MAG','REF','ref_new']
df_ctgs_melted = pd.melt(df_ctgs, id_vars=['Bin ID','ref_new'],
                    value_vars=['MAG', 'REF'],
                    var_name='Genome', value_name='Count')
df_ctgs_melted['category'] = df_ctgs_melted['Genome'] + '_' + df_ctgs_melted['ref_new']
df_ctgs_melted['log count'] = np.log10(df_ctgs_melted['Count'])

## tRNAs
df_trna = merged_hq_only[['Bin ID','Unique tRNAs_x','Unique tRNAs_y','ref_new']]
df_trna.columns = ['Bin ID','MAG','REF','ref_new']
df_trna_melted = pd.melt(df_trna, id_vars=['Bin ID', 'ref_new'],
                         value_vars=['MAG', 'REF'],
                         var_name='Genome', value_name='Count')
df_trna_melted['category'] = df_trna_melted['Genome'] + '_' + df_trna_melted['ref_new']

## Plot boxplots with three categories (RefSeq, GenBank, canine MAGs)
df_melted = df_ctgs_melted #df_trna_melted df_ctgs_melted
var = 'log count' #'log count' 'Count'

df_melted['category'] = df_melted['category'].str.replace('MAG_GenBank','Shanghai Dog MAG')
df_melted['category'] = df_melted['category'].str.replace('MAG_RefSeq','Shanghai Dog MAG')

## Plotting
width_mm = 45
height_mm = 46
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)
ax.clear()

categories = ['Shanghai Dog MAG', 'REF_RefSeq','REF_GenBank']
color_palette = ['#1b9e77', '#a6761d','#e6ab02']

sns.boxplot(data=df_melted,
            x='category', y=var,
            ax=ax,
            width=0.8,
            palette=color_palette,
            linewidth=1,
            fliersize=2)

ax.set_title('')
ax.set_ylabel('counts')
ax.set_xlabel('')
ax.set_xticklabels([])
ax.tick_params(axis='x', bottom=False)
sns.despine(fig, trim=False)

plt.tight_layout()
#plt.show()
plt.savefig("/data/Projects/ShanghaiDogs/intermediate-outputs/figures/sp_MAG-vs-ref_ctgs_boxplot.svg")

## 1) Assess normality - reshape to long format + Shapiro test
df_melted = df_ctgs_melted #df_trna_melted df_ctgs_melted
var = 'log count' #'Count' 'log count'

categories = df_melted["category"].unique()
normality_results = {}

for category in categories:
    stat, p = shapiro(df_melted[df_melted["category"] == category][var])
    normality_results[category] = p
    print(f"Shapiro-Wilk test for {category}: p-value = {p:.10f}")

## 2) Compute Wilcoxon test for SHD vs GenBank & SHD vs RefSeq
## one-to-one comparison of ctgs, tRNAs, and ribosomal genes

genbank_data['log count MAG'] = np.log10(genbank_data['Nr contigs'])
genbank_data['log count REF'] = np.log10(genbank_data['Number'])
refseq_data['log count MAG'] = np.log10(refseq_data['Nr contigs'])
refseq_data['log count REF'] = np.log10(refseq_data['Number'])

# Perform Wilcoxon test
p_values = []
results = []

stat, p_value = wilcoxon(genbank_data['Unique tRNAs_x'],genbank_data['Unique tRNAs_y'],alternative='greater')
results.append({'Comparison': 'tRNAs MAG vs REFs', 'Statistic': stat, 'p-value': p_value, 'Group': 'GenBank'})
p_values.append(p_value)
stat, p_value = wilcoxon(refseq_data['Unique tRNAs_x'], refseq_data['Unique tRNAs_y'],alternative='greater')
results.append({'Comparison': 'tRNAs MAG vs REFs', 'Statistic': stat, 'p-value': p_value, 'Group': 'RefSeq'})
p_values.append(p_value)
stat, p_value = wilcoxon(genbank_data['log count MAG'],genbank_data['log count REF'],alternative='less')
results.append({'Comparison': 'ctgs MAG vs REFs', 'Statistic': stat, 'p-value': p_value, 'Group': 'GenBank'})
p_values.append(p_value)
stat, p_value = wilcoxon(refseq_data['log count MAG'], refseq_data['log count REF'],alternative='less')
results.append({'Comparison': 'ctgs MAG vs REFs', 'Statistic': stat, 'p-value': p_value, 'Group': 'RefSeq'})
p_values.append(p_value)
stat, p_value = wilcoxon(genbank_data['16S rRNA_x'],genbank_data['16S rRNA_y'],alternative='greater')
results.append({'Comparison': '16S MAG vs REFs', 'Statistic': stat, 'p-value': p_value, 'Group': 'GenBank'})
p_values.append(p_value)
stat, p_value = wilcoxon(refseq_data['16S rRNA_x'], refseq_data['16S rRNA_y'],alternative='greater')
results.append({'Comparison': '16S MAG vs REFs', 'Statistic': stat, 'p-value': p_value, 'Group': 'RefSeq'})
p_values.append(p_value)
stat, p_value = wilcoxon(genbank_data['23S rRNA_x'],genbank_data['23S rRNA_y'],alternative='greater')
results.append({'Comparison': '23S MAG vs REFs', 'Statistic': stat, 'p-value': p_value, 'Group': 'GenBank'})
p_values.append(p_value)
stat, p_value = wilcoxon(refseq_data['23S rRNA_x'], refseq_data['23S rRNA_y'],alternative='greater')
results.append({'Comparison': '23S MAG vs REFs', 'Statistic': stat, 'p-value': p_value, 'Group': 'RefSeq'})
p_values.append(p_value)

# Correct for multiple testing (e.g., Benjamini-Hochberg)
_, p_values_corrected, _, _ = multipletests(p_values, method='fdr_bh')

# Add corrected p-values to the results DataFrame
for idx, result in enumerate(results):
    result['Corrected p-value'] = p_values_corrected[idx]

# Create DataFrame
comparison_results = pd.DataFrame(results)
comparison_results['Corrected p-value'] = comparison_results['Corrected p-value'].apply(lambda x: f"{x:.2e}")
