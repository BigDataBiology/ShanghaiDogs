import pandas as pd
from glob import glob

for v in ['KEGG_ko', 'CAZy', 'KEGG_Module', 'COG_category']:
    partials = glob(f'outputs/functional/*{v}.tsv.gz')
    partials.sort()
    data = pd.concat([pd.read_csv(f, sep='\t', index_col=0).squeeze() for f in partials], axis=1)
    data.index.name = v
    data.to_csv(f'outputs/functional_tables/{v}.tsv.gz', sep='\t', index=True, header=True)
    print(f'Saved {v} ...', data.shape)
