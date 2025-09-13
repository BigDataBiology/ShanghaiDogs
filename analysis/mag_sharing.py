import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

data = pd.read_csv('intermediate-outputs/05_dereplication/01_drep/ANI_9999/data_tables/Ndb.csv')
mimag_tab = pd.read_csv('data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv',index_col=0)
metadata = pd.read_csv('data/ShanghaiDogsMetadata/SH_Dog_metadata_red.csv', index_col=0)
metadata = metadata[~metadata.index.isin(['D000'])] # remove biological replicate

# Store households in a dictionary
household = metadata['Household.1'].to_dict()

# Get the original_id
data['ref_sample'] = data.reference.str.split('.').str[0].str.split('_').str[-1]
data['qry_sample'] = data.querry.str.split('.').str[0].str.split('_').str[-1]
data = data[~data['ref_sample'].isin(['D000'])] # remove biological replicate
data = data[~data['qry_sample'].isin(['D000'])] # remove biological replicate

# Filter MAG comparisons
data = data.query('alignment_coverage > 0.50')
data = data.query('ani > 0.95')

# Create sample pairs comparison dataframe
paired = data.groupby(['ref_sample', 'qry_sample'])['ani'].apply(list)
paired = paired.reset_index()
paired['shared_sp'] = paired['ani'].apply(lambda x: sum(x > 0.95 for x in x))
paired['shared_st'] = paired['ani'].apply(lambda x: sum(x > 0.99 for x in x))
paired = paired[['ref_sample', 'qry_sample', 'shared_sp', 'shared_st']]

paired.eval('fraction_shared_st = shared_st/shared_sp', inplace=True)
paired['is_same_household'] = paired[['ref_sample', 'qry_sample']].apply(lambda x: household[x[0]] == household[x[1]], axis=1)
paired = paired.query('ref_sample < qry_sample') # remove same animal comparison

### Create a subset (at least 10 shared species) for visualization
paired_10 = paired.query('shared_sp >= 10')
m_v = stats.mannwhitneyu(paired_10.query('is_same_household')['fraction_shared_st'], paired_10.query('~is_same_household')['fraction_shared_st'])

width_mm = 75
height_mm = 75
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)  # adjust size as needed
ax.clear()
sns.boxplot(data=paired_10, x='is_same_household', y='fraction_shared_st', ax=ax,
            boxprops={'facecolor':'None'}, showfliers=False, width=0.7)
sns.stripplot(data=paired_10, x='is_same_household', y='fraction_shared_st', ax=ax,
              alpha=0.5, size=2.5, jitter=0.2, palette=['#D3D3D3','#e7298a'])

ax.set_title('Bacterial strains', fontsize=11, pad=22)
ax.text(0.5, 1.1, f'p-value = {m_v.pvalue:.1e}', ha='center', va='center', transform=ax.transAxes)
ax.set_ylim(0,1)

ax.set_ylabel('Fraction shared strains', fontsize=10)
ax.set_xlabel('', fontsize=10)
ax.set_xticklabels(['Different\nHousehold', 'Same\nHousehold'],rotation=0,ha='center')

sns.despine()
plt.tight_layout()
plt.show()
fig.savefig('intermediate-outputs/figures/mag_sharing_ANI.svg')

print(f'''Mann-Whitney U test: p-value = {m_v.pvalue}''')

## Find bacterial species that are commonly >99% identical across SH dog cohort
## Find taxonomy of highly-shared MAGs
data['reference'] = data['reference'].str.replace('.fa.gz','')
data['querry'] = data['querry'].str.replace('.fa.gz','')
data_tax = pd.merge(data,mimag_tab[['Classification','Original ID']],left_on='reference',right_on='Original ID')
data_tax = pd.merge(data_tax,mimag_tab[['Classification','Original ID']],left_on='querry',right_on='Original ID')
data_tax.rename(columns={'Classification_x': 'Ref Classification', 'Classification_y': 'Qry Classification'}, inplace=True)
data_tax.drop(['Original ID_x','Original ID_y'],axis=1,inplace=True)

# Remove comparisons with eachself & for comparisons in the two directions, keep one
mask = data_tax['reference'] != data_tax['querry']
data_tax = data_tax[mask]
data_tax['Comparison'] = data_tax[['reference', 'querry']].apply(lambda x: tuple(sorted(x)), axis=1)
data_tax = data_tax.drop_duplicates('Comparison')
data_tax.drop(columns=['Comparison'], inplace=True)

# Proportion of same strain within the same species (same sp defined as ANI>95%)
data_tax_99 = data_tax[data_tax['ani']>0.99]
tax_grouping_99 = data_tax_99['Ref Classification'].value_counts().reset_index()
tax_grouping = data_tax['Ref Classification'].value_counts().reset_index()

merged_tax = pd.merge(tax_grouping_99,tax_grouping,right_on='Ref Classification',left_on='Ref Classification')
#merged_tax = merged_tax.query('count_y > 100')
merged_tax['comparisons >99 ANI'] = merged_tax['count_x']/merged_tax['count_y']*100
merged_tax.columns = ['Classification','ani >=99%','ani >=95%','% ani >=99% if same sp']

merged_tax[['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']] = \
    merged_tax['Classification'].str.split(';', expand=True)

merged_tax_sorted = merged_tax.sort_values(by='% ani >=99% if same sp', ascending=False)

# Create subsets
high_st_sharing = merged_tax_sorted[merged_tax_sorted['% ani >=99% if same sp']>70]
low_st_sharing = merged_tax_sorted[merged_tax_sorted['% ani >=99% if same sp']<10]

fig,ax = plt.subplots(figsize=(8,4))
order_colors = {phylum: plt.cm.Dark2(i) for i, phylum in enumerate(merged_tax_sorted['Phylum'].unique())}
ax.bar(high_st_sharing['Species'], high_st_sharing['% ani >=99% if same sp'], \
        color=[order_colors[order] for order in high_st_sharing['Phylum']])
ax.set_xticklabels(high_st_sharing['Species'],rotation=30, ha='right')
ax.set_ylim(0,100)
ax.set_ylabel('% comparisons with ANI >= 99%')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
fig.tight_layout()
#plt.show()
fig.savefig('figures/species_high_ANI.svg')

## Filter a specific genus only & plot boxplots
data_tax_genus = data_tax[data_tax['Ref Classification'].str.contains('Blautia')]
species_counts = data_tax_genus['Ref Classification'].value_counts()
print(species_counts)
species_to_keep = list(species_counts[species_counts >= 100].index)
species_to_keep.remove('d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Blautia_A;s__')
data_tax_genus_filt = data_tax_genus[data_tax_genus['Ref Classification'].isin(species_to_keep)]
data_tax_genus_filt = data_tax_genus_filt[~(data_tax_genus_filt['Ref Classification'].str.contains('Blautia_A'))]
data_tax_genus_filt.loc[:, 'Species'] = data_tax_genus_filt['Ref Classification'].str.split(';', expand=True)[6]
data_tax_genus_filt.loc[:, 'Species'] = data_tax_genus_filt['Species'].str.replace('s__','')
data_tax_genus_filt['is_same_household'] = data_tax_genus_filt[['ref_sample', 'qry_sample']].apply(lambda x: household[x[0]] == household[x[1]], axis=1)

data_tax_genus_filt['ani'] = data_tax_genus_filt['ani']*100

# Boxplot specific genus
fig, ax = plt.subplots()
ax.clear()

sns.boxplot(data=data_tax_genus_filt, x='Species', y='ani', ax=ax, boxprops={'facecolor':'None'}, showfliers=False)
sns.stripplot(data=data_tax_genus_filt[data_tax_genus_filt['is_same_household'] == False],
              x='Species', y='ani', ax=ax, alpha=0.3, color='darkgray')
sns.stripplot(data=data_tax_genus_filt[data_tax_genus_filt['is_same_household'] == True],
              x='Species', y='ani', ax=ax, alpha=0.6, color='green')

sns.despine()
ax.set_xticklabels(ax.get_xticklabels(), rotation=45,horizontalalignment='right')
ax.set_ylabel('ANI')
ax.set_xlabel('')
fig.tight_layout()
#plt.show()
fig.savefig('figures/Blautia_ANI.svg')

## Filter the most prevalent species & plot boxplots

## List of prevalent species 30 assembled MAGs / 50 dogs
prev_sp = ['Blautia hansenii', 'Faecalimonas umbilicata', 'Ruminococcus_B gnavus',
           'Enterocloster sp001517625', 'Clostridium_Q sp000435655', 'Collinsella intestinalis',
           'Blautia sp000432195', 'Megamonas funiformis','Amedibacillus dolichus',
           'Blautia_A sp900541345', 'Phocaeicola sp900546645', 'Schaedlerella glycyrrhizinilytica_A',
           'Fusobacterium_B sp900541465', 'Blautia_A caecimuris', 'Sutterella wadsworthensis_A',
           'Amedibacterium intestinale', 'Faecalimonas sp900550235', 'Phocaeicola coprocola',
           'Oliverpabstia sp000432335', 'Faecalibacterium sp900540455', 'Ventrimonas sp900538475',
           'Eisenbergiella sp900539715', 'Bacteroides sp900766005']

# Filtering
data_tax.loc[:, 'Species'] = data_tax['Ref Classification'].str.split(';', expand=True)[6]
data_tax['Species'] = data_tax['Species'].str.replace('s__','')
data_tax_filt = data_tax[data_tax['Species'].isin(prev_sp)]
comparisons_count = data_tax_filt['Ref Classification'].value_counts()
print(comparisons_count)

data_tax_filt['is_same_household'] = data_tax_filt[['ref_sample', 'qry_sample']].apply(lambda x: household[x[0]] == household[x[1]], axis=1)
data_tax_filt['ani'] = data_tax_filt['ani']*100

# Boxplots for prevalent species
width_mm = 180
height_mm = 95
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)
ax.clear()

sns.stripplot(data=data_tax_filt, hue='is_same_household',dodge=True,
              x='Species', y='ani', ax=ax, alpha=0.3, size=2, jitter=0.2, zorder=1)
sns.boxplot(data=data_tax_filt, x='Species', y='ani', ax=ax, boxprops={'facecolor':'None'},
            showfliers=False, dodge=True, hue='is_same_household',zorder=2)

# Remove duplicate legends and format plot
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:2], ['Same household', 'Different household'], title='')

sns.despine(fig, trim=False)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45,horizontalalignment='right')
ax.set_ylabel('ANI')
ax.set_xlabel('')
fig.tight_layout()

#plt.show()
fig.savefig('figures/prev_sp_MAG_sharing.svg')

### Plot a boxplot per each prevalent species
for species in data_tax_filt['Species'].unique():
    # Filter the data for the current species
    species_data = data_tax_filt[data_tax_filt['Species'] == species]

    # Calculate if significant
    m_v = stats.mannwhitneyu(species_data[species_data['is_same_household']]['ani'],
                             species_data[~species_data['is_same_household']]['ani'])
    print(f'{species} p-value: {m_v.pvalue}')

    # Create a new figure for each species
    width_mm = 70
    height_mm = 80
    figsize_inch = (width_mm / 25.4, height_mm / 25.4)
    fig, ax = plt.subplots(figsize=figsize_inch)
    sns.boxplot(data=species_data, x='is_same_household', y='ani', ax=ax,
                boxprops={'facecolor': 'None'}, showfliers=False)
    sns.stripplot(data=species_data, x='is_same_household', y='ani',
                  dodge=True, ax=ax, alpha=0.3, marker='o', size=2.5, jitter=0.2)

    # Add significance
    ax.text(0.5, 0.99, f'p-value = {m_v.pvalue:.1e}', ha='center', va='center',
            transform=ax.transAxes)
    # Labeling and styling
    ax.tick_params(labelsize=10)
    ax.set_title(f'{species}',fontsize=10,pad=10)
    ax.set_xlabel('')
    ax.set_ylabel('ANI (%)')
    ax.set_ylim([95,100])
    ax.set_xticklabels(['Between\nHouseholds','Within \nHouseholds'])
    sns.despine(fig, trim=False)

    fig.tight_layout()
    plt.show()
    #fig_name = 'figures/'+species.replace(' ', '_')+'__mag_sharing.svg'
    #fig.savefig(fig_name)
