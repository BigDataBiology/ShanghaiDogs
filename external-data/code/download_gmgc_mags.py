from jug import TaskGenerator
from os import path
import pandas as pd

@TaskGenerator
def download_gmgc_mags(species, habitat, country):
    import requests
    gmbc = pd.read_table('GMGC10.data/GMBC10.meta.tsv', index_col=0)
    sample_meta = pd.read_table('GMGC10.data/metadata/GMGC10.sample.meta.tsv.gz', index_col=1)
    gmbc = gmbc.query('category != "low-quality"')
    gmbc = gmbc[gmbc['GTDB_tk'].str.split(';').str[-1].str.contains(species)]
    gmbc['sample'] = gmbc['genome'].str.split('.').str[0]
    gmbc = gmbc.merge(sample_meta, left_on='sample', right_index=True)

    selected = gmbc
    if habitat is not None:
        selected = selected.query('habitat == @habitat')
    if country is not None:
        selected = selected.query('country == @country')
    for g in selected.index:
        ofile = f'GMGC_mags/{g}.fna.gz'
        if not path.exists(ofile):
            url = f'https://gmgc.embl.de/dumpFile.cgi?gnf={g}'
            r = requests.get(url, allow_redirects=True, stream=True)
            with open(ofile, 'wb') as f:
                while ch := r.raw.read(8192):
                    f.write(ch)
    return selected


download_gmgc_mags('Prevotella copri', 'dog gut', None)
download_gmgc_mags('Prevotella copri', 'cat gut', None)
download_gmgc_mags('Prevotella copri', 'human gut', 'Spain')
download_gmgc_mags('Prevotella copri', 'human gut', 'China')
download_gmgc_mags('Prevotella copri', 'human gut', 'United States of America')
