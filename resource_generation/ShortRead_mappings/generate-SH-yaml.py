import yaml
import glob
from collections import defaultdict

BASEDIR = '../../data/ShanghaiDogsFastQ/ILM/'

samples = glob.glob(f'{BASEDIR}/*')

r = {}
for s in samples:
    fqs = glob.glob(f'{s}/*')
    fqs.sort()
    fqs = [f.removeprefix(BASEDIR) for f in fqs]
    [f1,f2] = fqs
    sample = s.removeprefix(BASEDIR)
    r[sample] = [{'paired': [f1,f2]}]

r = {
        'basedir': BASEDIR,
        'samples' : dict(r),
    }

with open('SH_dogs.yaml', 'wt') as out:
    out.write(yaml.dump(r))
