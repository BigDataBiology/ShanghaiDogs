import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pysankey2 import Sankey

os.chdir('/data/Projects/ShanghaiDogs')
plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

# Import files
polish_eval = pd.read_csv("intermediate-outputs/polishing_evaluation/results/Checkm2_final_polish_eval.csv", \
                  delimiter=',',header=0,index_col=0)

# Reformat 'Cov' categories
polish_eval.loc[polish_eval['Mean_cov LR'] < 40, 'Cov'] = 'LOW'
polish_eval.loc[polish_eval['Mean_cov LR'] > 40, 'Cov'] = 'HIGH'

# Order by quality evolution through polishing
order = ['High-to-High','Medium-to-High','Medium-to-Medium','Low-to-Medium','High-to-Medium']

polish_eval['Qual-evol'] = pd.Categorical(polish_eval['Qual-evol'],
                                          categories=order, ordered=True)
polish_eval['Qual-evol'].dropna(inplace=True)
polish_eval = polish_eval.sort_values(by='Qual-evol')

## Scatterplot of diff in completeness (Polca - Flye)
# PLOTTING RESULTS
width_mm = 160 # 50 and 75
height_mm = 80
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig, ax = plt.subplots(figsize=figsize_inch)
ax.axhline(0, color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(
    x=polish_eval['Mean_cov LR'],
    y=polish_eval['Diff_Flye-Polca'],
    palette=['#11654c','#1b9e77','#e6ab02','#8f6b03','#d95f02'],
    hue=polish_eval['Qual-evol'],
    s=10,
    alpha=1,
    zorder=2)

ax.set_xlim(0,400)
ax.set_ylabel('Completeness variation (%)')
ax.set_xlabel('MAG coverage')
ax.legend_.set_title("")
sns.despine()
plt.tight_layout()

#plt.show()
plt.savefig('intermediate-outputs/figures/scatterplot_polish_eval_0-400_cov.svg')

## Boxplot of completeness diff at each polishing step (low vs high)
polish_eval_low = polish_eval[polish_eval['Cov'] == 'LOW']
polish_eval_low = polish_eval_low.reset_index()
polish_eval_high = polish_eval[polish_eval['Cov'] == 'HIGH']
polish_eval_high = polish_eval_high.reset_index()
print(polish_eval_low['Cov'])
print(polish_eval_high['Cov'])

order = ['Diff_Flye-Med','Diff_Flye-Poly','Diff_Flye-Polca']
polish_cov_long = pd.melt(polish_eval_high, id_vars=['index'],
                    value_vars=['Diff_Flye-Med', 'Diff_Flye-Poly','Diff_Flye-Polca'],
                    var_name='Polishing', value_name='Difference')
polish_cov_long['Polishing'] = pd.Categorical(polish_cov_long['Polishing'],
                                                   categories=order, ordered=True)
polish_cov_long = polish_cov_long.sort_values(by='Polishing')

# PLOTTING RESULTS
width_mm = 50
height_mm = 50
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

# Striplot

fig, ax = plt.subplots(figsize=figsize_inch)

# Adding individual data points using stripplot with increased jitter
sns.stripplot(data=polish_cov_long,
              x='Polishing', y='Difference',
              ax=ax,
              hue='Polishing',
              palette="YlOrBr",
              size=1,
              alpha=0.5,  # Transparency of the dots
              jitter=0.3,  # Increase jitter to spread out the dots
              order=order)

sns.boxplot(data=polish_cov_long,
            x='Polishing', y='Difference',
            ax=ax,
            width=0.8,
            palette=[(1, 1, 1, 0) for _ in "YlOrBr"], # Transparent boxplot
            linewidth=1,
            fliersize=0,  # Hide the outlier markers
            order=order)

ax.legend_.remove()
ax.set_xlabel('')
ax.set_xticklabels([])
ax.set_ylim(-20, 20)
ax.tick_params(axis='x', bottom=False)
sns.despine(fig, trim=False)

plt.tight_layout()
#plt.show()
plt.savefig('intermediate-outputs/figures/boxplot_compl_HIGH_cov.svg')

## Boxplot of contamination diff at each polishing step (low vs high)
polish_eval_low = polish_eval[polish_eval['Cov'] == 'LOW']
polish_eval_low = polish_eval_low.reset_index()
polish_eval_high = polish_eval[polish_eval['Cov'] == 'HIGH']
polish_eval_high = polish_eval_high.reset_index()

print(polish_eval_low['Cov'])
print(polish_eval_high['Cov'])

order = ['Contam_diff_Flye_Med','Contam_diff_Flye_Poly','Contam_diff_Flye_Polca']


polish_cov_long = pd.melt(polish_eval_low, id_vars=['index'],
                    value_vars=['Contam_diff_Flye_Med', 'Contam_diff_Flye_Poly','Contam_diff_Flye_Polca'],
                    var_name='Polishing', value_name='Difference')
polish_cov_long['Polishing'] = pd.Categorical(polish_cov_long['Polishing'],
                                                   categories=order, ordered=True)
polish_cov_long = polish_cov_long.sort_values(by='Polishing')

# PLOTTING RESULTS
width_mm = 50
height_mm = 50
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

# Striplot

fig, ax = plt.subplots(figsize=figsize_inch)

# Adding individual data points using stripplot with increased jitter
sns.stripplot(data=polish_cov_long,
              x='Polishing', y='Difference',
              ax=ax,
              hue='Polishing',
              palette="YlOrBr",
              size=1,
              alpha=0.5,  # Transparency of the dots
              jitter=0.3,  # Increase jitter to spread out the dots
              order=order)

sns.boxplot(data=polish_cov_long,
            x='Polishing', y='Difference',
            ax=ax,
            width=0.8,
            palette=[(1, 1, 1, 0) for _ in "YlOrBr"], # Transparent boxplot
            linewidth=1,
            fliersize=0,  # Hide the outlier markers
            order=order)

ax.legend_.remove()
ax.set_xlabel('')
ax.set_xticklabels([])
ax.set_ylim(-10, 10)
ax.tick_params(axis='x', bottom=False)
sns.despine(fig, trim=False)

plt.tight_layout()
#plt.show()
plt.savefig('intermediate-outputs/figures/boxplot_contam_LOW_cov.svg')

# Create Sankey diagrams for Quality evolution
# For MAGs with high coverage (>40X)
polish_eval_low = polish_eval[polish_eval['Cov']=='LOW'] # <40X cov n=1205
polish_eval_high = polish_eval[polish_eval['Cov']=='HIGH'] # >40X cov n=1361

# Plot Sankey diagrams
df = polish_eval_low.loc[:,['Flye_qual','Medaka_qual','Poly_qual','Polca_qual']]
df.columns = ['layer1','layer2','layer3','layer4']

layer_labels= {'layer1':['high_quality','medium_quality','low_quality'],
               'layer2':['high_quality','medium_quality','low_quality'],
               'layer3':['high_quality','medium_quality','low_quality'],
               'layer4':['high_quality','medium_quality','low_quality']}
global_color_palette = {'high_quality': '#1b9e77','medium_quality': '#e6ab02','low_quality':'#d95f02'}

sky_auto_global_colors = (
    Sankey(df,
           colorMode="global",
           layerLabels = layer_labels,
           colorDict=global_color_palette,
           stripColor='left'))

width_mm = 60
height_mm = 50
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

fig,ax = sky_auto_global_colors.plot(figSize=figsize_inch,
                                     fontSize=10)
#plt.tight_layout()
plt.savefig('intermediate-outputs/figures/sankey_low_cov.svg')
plt.show()