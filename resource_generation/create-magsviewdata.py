from os import makedirs
import pandas as pd
import json

data = pd.read_csv('../intermediate-outputs/07_ribosomal_genes/microbe-atlas.csv.gz')
data = data.groupby('MAG').apply(lambda x: x.T.to_dict()).to_dict()
magdata = {}
for mag,vals in data.items():
    mag16s = [c for _,c in vals.items()]
    for m in mag16s:
        del m['MAG']
    magdata[mag] = {'16S': mag16s}

makedirs('../intermediate-outputs/magsviews/genome-data', exist_ok=True)
for key,mdata in magdata.items():
    with open(f'../intermediate-outputs/magsviews/genome-data/{key}.json', 'w') as f:
        json.dump(mdata, f)
