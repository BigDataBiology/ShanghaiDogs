import pandas as pd
from glob import glob

kegg_kos = glob('outputs/functional/*KEGG_ko.tsv.gz')
data = pd.concat([pd.read_csv(f, sep='\t', index_col=0).squeeze() for f in kegg_kos], axis=1)
data.index.name = 'KO'
data.to_csv('outputs/KEGG_ko.tsv.gz', sep='\t', index=True, header=True)
