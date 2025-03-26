import yaml
import glob
from collections import defaultdict
from jug import Task



@Task
def generate_yarlagadda_yaml():
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


@Task
def generate_berlin_dogs_yaml():
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

@Task
def generate_shanghai_dogs_yaml():
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

@Task
def generate_USA_pets_Nomnomnow():
    BASEDIR = '../../external-data/data/USA_pets_Nomnomnow/'
    fqs= glob.glob(f'{BASEDIR}*.gz')
    r = {}
    for f in fqs:
        f = f.removeprefix(BASEDIR)
        sample = f.removesuffix('.fastq.gz')
        r[sample] = [{'single': [f]}]
    r = {
            'basedir': BASEDIR,
            'samples': r,
    }
    print(r)
    with open('USA_pets_Nomnomnow.yaml', 'wt') as out:
        out.write(yaml.dump(r))

@Task
def generate_Coelho_2018():
    BASEDIR = '../../external-data/data/Coelho_dogs_2018/'

    samples = glob.glob(f'{BASEDIR}/*')

    r = {}
    for s in samples:
        fqs = glob.glob(f'{s}/*')
        fqs.sort()
        f1s = [f for f in fqs if f.endswith('.1.fq.gz')]
        assert len(f1s) == len(fqs) // 2
        f1s = [f.removeprefix(BASEDIR) for f in f1s]
        reads = []
        for f1 in f1s:
            f2 = f1.replace('.1.fq.gz', '.2.fq.gz')
            reads.append({'paired': [f1,f2]})
        sample = s.removeprefix(BASEDIR)
        if reads:
            r[sample] = reads

    r = {
            'basedir': BASEDIR,
            'samples': dict(r),
        }

    with open('Coelho_2018.yaml', 'wt') as out:
        out.write(yaml.dump(r))


@Task
def generate_Allaway():
    BASEDIR = '../../external-data/data/Allaway/WGS_RANDOM/'

    samples = glob.glob(f'{BASEDIR}/*')

    r = {}
    for s in samples:
        fqs = glob.glob(f'{s}/*')
        fqs.sort()
        f1s = [f for f in fqs if f.endswith('.1.fq.gz')]
        assert len(f1s) == len(fqs) // 2
        f1s = [f.removeprefix(BASEDIR) for f in f1s]
        reads = []
        for f1 in f1s:
            f2 = f1.replace('.1.fq.gz', '.2.fq.gz')
            reads.append({'paired': [f1,f2]})
        sample = s.removeprefix(BASEDIR)
        if reads:
            r[sample] = reads

    r = {
            'basedir': BASEDIR,
            'samples': dict(r),
        }

    with open('Allaway.yaml', 'wt') as out:
        out.write(yaml.dump(r))
