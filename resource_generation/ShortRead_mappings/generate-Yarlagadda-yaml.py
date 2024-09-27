import yaml
import glob
from collections import defaultdict

BASEDIR = '../../external-data/data/Yarlagadda/'

fs = glob.glob(f'{BASEDIR}/*/*.gz')
f1s = [f for f in fs if '.1.fq.gz' in f]
r = defaultdict(list)
for f1 in f1s:
    f1 = f1.removeprefix(BASEDIR)
    f2 = f1.replace('.1.fq.gz', '.2.fq.gz')
    sample = f1.split('/')[0]
    r[sample].append({'paired': [f1,f2]})

r = {
        'basedir': BASEDIR,
        'samples' : dict(r),
    }

with open('Yarlagadda_dogs.yaml', 'wt') as out:
    out.write(yaml.dump(r))
