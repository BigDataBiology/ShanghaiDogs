import pandas as pd
from glob import glob
import yaml

shd_mapstats = glob('outputs/*_shd_mapstats.txt')
shd_sp_mapstats = glob('outputs/*_shd_sp_mapstats.txt')
shd_nc_mapstats = glob('outputs/*_shd_nc_mapstats.txt')
shd_sp_nc_mapstats = glob('outputs/*_shd_sp_nc_mapstats.txt')

def load_all(files):
    all_ = {}
    for file in sorted(files):
        sample = file.split('_')[0].split('/')[-1]
        all_[sample] = pd.read_csv(file, sep='\t', index_col=0).squeeze()[['total', 'aligned']]
    return pd.DataFrame.from_dict(all_, orient='index')

aligned = load_all(shd_mapstats)
aligned_sp = load_all(shd_sp_mapstats)
aligned_nc = load_all(shd_nc_mapstats)
aligned_sp_nc = load_all(shd_sp_nc_mapstats)


data = pd.DataFrame(
               {'total': aligned['total']
                ,'aligned': aligned['aligned']
                ,'total_sp': aligned_sp['total']
                ,'aligned_sp': aligned_sp['aligned']
                ,'total_nc': aligned_nc['total']
                ,'aligned_nc': aligned_nc['aligned']
                ,'total_sp_nc': aligned_sp_nc['total']
                ,'aligned_sp_nc': aligned_sp_nc['aligned']
              }) \
    .reset_index() \
    .rename(columns={'index': 'sample'})

assert data.dropna().eval('total == total_sp').all()
data.drop('total_sp', axis=1, inplace=True)

assert data.dropna().eval('total == total_nc').all()
data.drop('total_nc', axis=1, inplace=True)

assert data.dropna().eval('total == total_sp_nc').all()
data.drop('total_sp_nc', axis=1, inplace=True)

nestle = set(line.strip() for line in open('../../external-data/data/Coelho_dogs_2018/PRJEB20308'))
berlin = set(s.split('_')[0] for s in yaml.safe_load(open('berlin_dogs.yaml'))['samples'])
Yarlagadda = set(yaml.safe_load(open('Yarlagadda_dogs.yaml'))['samples'])
shanghai = set(yaml.safe_load(open('SH_dogs.yaml'))['samples'])
nomnomnow = set(yaml.safe_load(open('USA_pets_Nomnomnow.yaml'))['samples'])
allaway = set(yaml.safe_load(open('Allaway.yaml'))['samples'])

def group_for(s):
    if s in nestle:
        return 'Nestlé'
    if s in berlin:
        return 'Berlin'
    if s in Yarlagadda:
        return 'Yarlagadda'
    if s in shanghai:
        return 'Shanghai'
    if s in nomnomnow:
        return 'NomNomNow'
    if s in allaway:
        return 'Allaway'
data['group'] = data['sample'].map(group_for)
data.to_csv('../../intermediate-outputs/external_datasets_mappings/reads_mapped_shd.tsv', sep='\t', index=False)
