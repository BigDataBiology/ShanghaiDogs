import yaml
import glob
from collections import defaultdict

BASEDIR = '../../external-data/data/berlin_dogs/'

fs = glob.glob(f'{BASEDIR}*.gz')
f1s = [f for f in fs if '_1.fq.gz' in f]
r = defaultdict(list)
for f1 in f1s:
    f1 = f1.removeprefix(BASEDIR)
    f2 = f1.replace('_1.fq.gz', '_2.fq.gz')
    #animal_code = f1[:-len('_Lx_x.fq.gz')]
    animal_code = f1.split('-')[0]
    r[animal_code].append({'paired': [f1,f2]})

r = {
        'basedir': BASEDIR,
        'samples' : dict(r),
    }

with open('berlin_dogs.yaml', 'wt') as out:
    out.write(yaml.dump(r))
