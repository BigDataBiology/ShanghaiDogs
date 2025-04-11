import pandas as pd
from glob import glob
BASEDIR = '../intermediate-outputs/Prodigal/'
origin = pd.read_csv(f'{BASEDIR}/SHD.ORF.orig.tsv.xz', sep='\t', names=['SHD_ORF', 'Sample', 'ORF'], index_col=0)
origin['sample_orf'] = origin["Sample"] + '__' + origin["ORF"]
origin = origin.drop(['Sample', 'ORF'], axis=1)
origin = origin.reset_index().set_index('sample_orf')

clusters = pd.read_csv(f'{BASEDIR}/SHD.clusters.tsv.xz', sep='\t', names=['SHD_ORF', 'rel1', 'SHD.100NT', 'rel2', 'SHD.95NT'])
eqs = clusters.query('rel1 == "=" & rel2 == "="')

total_95nt = len(set(clusters['SHD.95NT']))
eq_95nt = len(set(eqs['SHD.95NT']))
total_95nt == eq_95nt
assert total_95nt == eq_95nt

eqs = eqs[['SHD.95NT', 'SHD_ORF']].sort_values(by=['SHD.95NT','SHD_ORF'])
eqs = eqs.drop_duplicates(subset=['SHD.95NT'], keep='first')
eqs.set_index("SHD_ORF", inplace=True)
origin = origin.join(eqs, on='SHD_ORF').dropna()
origin = origin[['SHD.95NT']]


emappers = sorted(glob('../intermediate-outputs/eggNOG_annot_contigs/D*/*.annotations'))
annotations = []
for e in emappers:
    sample = e.split('/')[3]
    c = pd.read_csv(e, sep='\t', skiprows=4, skipfooter=3, index_col=0, engine='python')
    c.rename(index= lambda i: sample + '__' + i, inplace=True)
    annotations.append(c)
annotations = pd.concat(annotations)
annotations = annotations.join(origin, how='inner')
annotations = annotations.reset_index().set_index("SHD.95NT").copy()
annotations.sort_index(inplace=True)
annotations.to_csv(f'{BASEDIR}/SHD.95NT.emapper.annotations.gz', sep='\t')
