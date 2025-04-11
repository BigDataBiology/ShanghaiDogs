# -*- coding: utf-8 -*-

"""
Created on Thu Feb 28 14:58:12 2024
@author: Anna Cusco
"""

import os
import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

os.chdir('/data/Projects/ShanghaiDogs')
plt.rcParams['svg.fonttype'] = 'none' #to avoid transforming the font to plot

# Import files
QC_all = "intermediate-outputs/00_quality_control/QC_ONT_raw_reads.txt"

with open(QC_all, 'r') as file:
  QC_all = pd.read_csv(file, delimiter='\t',header=None,index_col=0)

QC_all.columns = ['Total gigabases','N50 length','Median length','Median Q-score']
QC_all = QC_all.applymap(lambda x: re.sub(r'^.*?: ', '', x) if isinstance(x, str) else x)
QC_all = QC_all.applymap(pd.to_numeric)
QC_all.to_csv('data/ShanghaiDogsTables/ONT_raw_reads_report.csv')

### Basic histograms
width_mm = 60
height_mm = 55
figsize_inch = (width_mm / 25.4, height_mm / 25.4)

# Median Q-score
fig,ax = plt.subplots(figsize=figsize_inch)

QC_all.hist(column='Median Q-score', bins=10, grid=False,
            color='#a6761d', rwidth=0.85, ax=ax, alpha=0.7)

plt.title('')
plt.xlabel('Q score (median)',fontsize=9)
plt.ylabel('Number of samples',fontsize=9)
plt.xticks(fontsize=9)  # Set x-axis tick font size
plt.yticks(fontsize=9)  # Set y-axis tick font size

plt.tight_layout()
plt.show()
#plt.savefig("intermediate-outputs/figures/hist_median_quality.svg")

# Median Length
fig,ax = plt.subplots(figsize=figsize_inch)

QC_all.hist(column='Median length', bins=10, grid=False,
            color='#e6ab02', rwidth=0.85, ax=ax, alpha=0.7)

ax.set_xticks(np.arange(0, 14000, 4000))
ax.set_yticks(np.arange(0, 12, 2))

plt.title('')
plt.xlabel('Read length (median)',fontsize=9)
plt.ylabel('Number of samples',fontsize=9)
plt.xticks(fontsize=9)  # Set x-axis tick font size
plt.yticks(fontsize=9)  # Set y-axis tick font size

plt.tight_layout()
plt.show()
#plt.savefig("intermediate-outputs/figures/hist_median_length.svg")

# Total Gbps
fig,ax = plt.subplots(figsize=figsize_inch)

QC_all.hist(column='Total gigabases', bins=10, grid=False,
            color='#66a61e', rwidth=0.85, ax=ax, alpha=0.7)

ax.set_xticks(np.arange(20, 70, 10))

plt.title('')
plt.xlabel('Total Gbp',fontsize=9)
plt.ylabel('Number of samples',fontsize=9)
plt.xticks(fontsize=9)  # Set x-axis tick font size
plt.yticks(fontsize=9)  # Set y-axis tick font size

plt.tight_layout()
plt.show()
#plt.savefig("intermediate-outputs/figures/hist_total_Gbp.svg")
