from os import path
import gzip
import pandas as pd
from glob import glob
import os
TARGETDIR = '/home/luispedro/ShanghaiDogsMAGs'

def copy_rename(fafile, ix, extra):
    oname = f'{TARGETDIR}/SHD1_{ix:04}.fna.gz'
    contig_ix = 1
    with gzip.open(oname, 'wt') as ofile, \
           gzip.open(fafile, 'rt') as ifile:
       for line in ifile:
           if line[0] == '>':
               orig = line[1:].strip()
               line = f'>SHD1_{ix:04}_{contig_ix} {orig} {extra}\n';
               contig_ix +=1
           ofile.write(line)
    return oname

das_out = glob('*_DASTool_bins/*.fa.gz')
extra = glob('*_DASTool_bins/*subsets/*.fa.gz')
das_out.extend(extra)
das_out.sort()

allq = []
for q in glob('*_DASTool_bins/*_qual_report.csv'):
    q = pd.read_csv(q)
    if q.iloc[0]['sample_id'] != 'D024':
        q.rename(columns={'index':'Name'}, inplace=True)
    allq.append(q)
allq = pd.concat(allq)
allq.eval('Qscore = Completeness - 5* Contamination', inplace=True)
allq.sort_values(by='Qscore', inplace=True)
allq.eval('Qscore_neg = -Qscore', inplace=True)
allq.sort_values(by=['Qscore_neg', 'sample_id', 'Name'], inplace=True, )
allq = allq.query('Quality != "low_quality"')
allq = allq.rename(columns={'1': 'Nr_contigs'}, )


newnames = []
for _, r in allq.iterrows():
    fafile = f'{r["sample_id"]}_DASTool_bins/{r["Name"]}_{r["sample_id"]}.fa.gz'
    if str(r['Name']).startswith('10G_20G'):
        fafile = f'{r["sample_id"]}_DASTool_bins/subsets/{r["Name"]}_{r["sample_id"]}.fa.gz'
    if not path.exists(fafile):
        raise KeyError(f'Not found: {fafile}')
    cur = copy_rename(fafile, len(newnames), f'{r["Name"]} {r["sample_id"]}')
    newnames.append(cur)

allq['filename'] = newnames
columns_to_keep = ['filename', 'sample_id', 'Name', 'Completeness', 'Contamination', 'Quality']
allq = allq[columns_to_keep].rename(columns={'filename':'Filename', 'sample_id': 'Sample'})
allq.to_csv(f'{TARGETDIR}/mag_meta.tsv.gz', sep='\t', index=False)
