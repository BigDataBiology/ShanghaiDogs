from os import makedirs
import pandas as pd
import json
import argnorm.lib

ribo_data = pd.read_csv('../intermediate-outputs/07_ribosomal_genes/microbe-atlas.csv.gz')
ribo_data = ribo_data.groupby('MAG').apply(lambda x: x.T.to_dict()).to_dict()
magdata = {}
for mag,vals in ribo_data.items():
    mag16s = [c for _,c in vals.items()]
    for m in mag16s:
        del m['MAG']
    magdata[mag] = {'16S': mag16s}

args_mags = pd.read_csv("../intermediate-outputs/06_ARG/MAGs-ARGs_ALL_filt.txt", index_col=0)
args_mags['Bin ID'] = args_mags["Bin ID"].str.split('.').str[0]
args_mags['ARO'] = args_mags['ARO'].map(lambda ar: f'ARO:{ar}')

resf = argnorm.lib.get_aro_mapping_table('resfinder')
in_resfinder = set(resf['ARO'].dropna().unique())
args_mags['InResfinder'] = args_mags['ARO'].map(in_resfinder.__contains__)

args_mags = args_mags[['Bin ID', 'Predicted_Protein', 'Best_Hit_ARO', 'Cut_Off', 'Best_Identities', 'Percentage Length of Reference Sequence', 'Drug Class', 'ARO', 'InResfinder']].rename(columns={
    'Predicted_Protein': 'Sequence',
    'Best_Hit_ARO': 'ArgName',
    'Best_Identities': 'Identity',
    'Percentage Length of Reference Sequence': 'Coverage',
})
args_mags = args_mags.groupby('Bin ID').apply(lambda g: g.drop('Bin ID', axis=1).to_dict(orient='records')).to_dict()

for key,vals in args_mags.items():
    # We want Perfect hits first, then sorted by Identity and Coverage
    vals.sort(key=lambda x: ((0 if x['Cut_Off'] == "Perfect" else -1), x['Identity'], x['Coverage']), reverse=True)
    magdata[key]['ARGs'] = vals

for vals in magdata.values():
    if 'ARGs' not in vals:
        vals['ARGs'] = []

makedirs('../intermediate-outputs/magsviews/genome-data', exist_ok=True)
for key,mdata in magdata.items():
    with open(f'../intermediate-outputs/magsviews/genome-data/{key}.json', 'w') as f:
        json.dump(mdata, f)
