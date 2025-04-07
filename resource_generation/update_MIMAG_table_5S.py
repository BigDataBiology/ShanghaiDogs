import sys
import pandas as pd
import os
import gzip

os.chdir('/data/Projects/ShanghaiDogs/')

# Import quality report & RNAmmer output file
SHD_qual = pd.read_csv('intermediate-outputs/tables/SHD_bins_MIMAG_report_old.csv',
                       sep=',', index_col=0)
rrnammer_out = pd.read_csv('intermediate-outputs/rnammer/rnammer_output.txt',
                           sep='\t', comment='#',header=None)

# Reformat output
SHD_qual.index = SHD_qual.index.str.replace('.fna.gz','')
rrnammer_out.columns = ['seqname','source','feature','start','end','score','+/-','frame','attribute','empty']
rrnammer_out.drop(['empty'],axis=1, inplace=True)
rrnammer_out['MAG ID'] = rrnammer_out['seqname'].str.replace(r'_\d+$', '', regex=True)

# Count rRNAs within each MAG
rrna_count = rrnammer_out.groupby(['MAG ID','attribute']).size().reset_index(name='count')

# From long to wide format
rrna_count_wide = rrna_count.pivot(index='MAG ID', columns='attribute', values='count').fillna(0)

# Link barrnap results with rrnammer results
predicted_ribosomals = pd.merge(rrna_count_wide,SHD_qual[['16S rRNA','23S rRNA','5S rRNA',
                                                          '16S partial','23S partial','5S partial']],
                                right_index=True,left_index=True)

SHD_qual_updated = pd.merge(SHD_qual,rrna_count_wide,
                            right_index=True,left_index=True,how='left')
SHD_qual_updated.fillna(0, inplace=True)

# If 5S is 0, update it with rrnammer results
for index in SHD_qual_updated.index:
    if SHD_qual_updated.loc[index, '5S rRNA'] == 0:
        SHD_qual_updated.loc[index, '5S rRNA'] = SHD_qual_updated.loc[index, '5s_rRNA']

mimag_counts_old = SHD_qual_updated['MIMAG'].value_counts()

SHD_qual_updated.drop(['5S partial','16s_rRNA','23s_rRNA','5s_rRNA','MIMAG'],inplace=True,axis=1)

# Updated MIMAG counts (ribosomals>0, and tRNAs>=18)
for index in SHD_qual_updated.index:
    if SHD_qual_updated.loc[index, 'Quality'] == 'high-quality':
        has_ribosomals = (
            SHD_qual_updated.loc[index, '16S rRNA'] >= 1 and
            SHD_qual_updated.loc[index, '23S rRNA'] >= 1 and
            SHD_qual_updated.loc[index, '5S rRNA'] >= 1
        )
        has_tRNAs = SHD_qual_updated.loc[index, 'Unique tRNAs'] >= 18

        if has_ribosomals and has_tRNAs:
            SHD_qual_updated.loc[index, 'MIMAG'] = 'Yes'
        else:
            SHD_qual_updated.loc[index, 'MIMAG'] = 'No'
    else:
        SHD_qual_updated.loc[index, 'MIMAG'] = 'No'

mimag_counts_new = SHD_qual_updated['MIMAG'].value_counts()
SHD_qual_updated.to_csv('data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv',sep=',')
