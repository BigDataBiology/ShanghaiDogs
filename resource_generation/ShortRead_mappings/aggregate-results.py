import pandas as pd
from glob import glob
import yaml

shd_mapstats = glob('outputs/*_shd_mapstats.txt')
shd_sp_mapstats = glob('outputs/*_shd_sp_mapstats.txt')

aligned = {}
for file in sorted(shd_mapstats):
    sample = file.split('_')[0].split('/')[-1]
    aligned[sample] = pd.read_csv(file, sep='\t', index_col=0).squeeze()[['total', 'aligned']]
aligned = pd.DataFrame.from_dict(aligned, orient='index')

aligned_sp = {}
for file in sorted(shd_sp_mapstats):
    sample = file.split('_')[0].split('/')[-1]
    aligned_sp[sample] = pd.read_csv(file, sep='\t', index_col=0).squeeze()[['total', 'aligned']]
aligned_sp = pd.DataFrame.from_dict(aligned_sp, orient='index')


data = pd.DataFrame(
               {'total': aligned['total']
                ,'aligned': aligned['aligned']
                ,'total_sp': aligned_sp['total']
                ,'aligned_sp': aligned_sp['aligned']
              }) \
    .reset_index() \
    .rename(columns={'index': 'sample'})

assert data.dropna().eval('total == total_sp').all()

data.drop('total_sp', axis=1, inplace=True)

nestle = set(line.strip() for line in open('../../external-data/data/Coelho_dogs_2018/PRJEB20308'))
berlin = set(s.split('_')[0] for s in yaml.safe_load(open('berlin_dogs.yaml'))['samples'])
Yarlagadda = set(yaml.safe_load(open('Yarlagadda_dogs.yaml'))['samples'])
shanghai = set(yaml.safe_load(open('SH_dogs.yaml'))['samples'])
nomnomnow = set(yaml.safe_load(open('USA_pets_Nomnomnow.yaml'))['samples'])
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
data['group'] = data['sample'].map(group_for)
data.to_csv('../../intermediate-outputs/external_datasets_mappings/reads_mapped_shd.tsv', sep='\t', index=False)
