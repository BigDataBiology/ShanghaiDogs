from collections import defaultdict
from glob import glob
import os 
from os import path
import pathlib

BATCH1_ILM = 'Public/ANNA/wangxinyu0616030/Release-X101SC23043778-Z01-J001-B1-6_20230615/NGSdata-X101SC23043778-Z01-J001/01.CleanData'
BATCH2_ILM = 'Public/ANNA/wangxinyu0616030/Release-X101SC23043778-Z01-J002-B1-6_20230619/NGSdata-X101SC23043778-Z01-J002/01.CleanData'
BATCH3_ILM = 'Public/Result-X101SC23043778-Z01-J005-B1-6_delivery_20230818/NGS/01.CleanData'
BATCH4_ILM = 'Public/X101SC23043778-Z01/Clean_Illumina_AUG_2023/X101SC23043778-Z01-F008'

found_data = defaultdict(list)

for base in [BATCH3_ILM, BATCH4_ILM]:
    for fq1 in glob(f'{base}/D0*/D0*_350.fq1.gz'):
        fq2 = fq1.replace('.fq1.', '.fq2.')
        if not path.exists(fq2):
            raise IOError(f'Expected {fq2}')
        sample = fq1.split('/')[-2]
        tech = 'ILM'
        assert sample.startswith('D0')
        found_data[sample,tech].append((fq1, fq2))
#        print(f'{sample} {tech}')
#        print(f'  {fq1}')
#        print(f'  {fq2}')
#        print()



BATCH1_ONT = 'Public/ANNA/wangxinyu0616030/Release-X101SC23043778-Z01-J001-B1-6_20230615/ONTdata_X101SC23043778-Z01-J001/'
BATCH2_ONT = 'Public/ANNA/wangxinyu0616030/Release-X101SC23043778-Z01-J002-B1-6_20230619/ONTdata_X101SC23043778-Z01-J002'
for base in [BATCH1_ONT, BATCH2_ONT]:
    for fq in glob(f'{base}/D0*/*/*/D0*fastq_pass.gz'):
        sample = fq.split('/')[-4]
        flowcell = fq.split('/')[-2]
        tech = 'ONT'
        assert sample.startswith('D0')
        if sample == 'D007':
            # this batch was low quality
            continue
        found_data[sample, tech].append(fq)
        #print(f'{sample} {tech}')
        #print(f'  {fq} {flowcell}')
        #print()

BATCH3_ONT = 'Public/Result-X101SC23043778-Z01-J005-B1-6_delivery_20230818/ONT'
for fq in glob(f'{BATCH3_ONT}/D0*/D0*fastq_pass.gz'):
    sample = fq.split('/')[-2]
    assert sample.startswith('D0')
    found_data[sample, 'ONT'].append(fq)
#    print(f'{sample} {tech}')
#    print(f'  {fq}')
#    print()

D007_SPECIAL_CONCAT = 'Public/X101SC23043778-Z01/D007_new/D007_pass.fastq.gz'

found_data['D007', 'ONT'] = [D007_SPECIAL_CONCAT]

for (s,t),fs in found_data.items():
    targetdir = f'../data/ShanghaiDogsFastQ/{t}/{s}'
    os.makedirs(targetdir, exist_ok=True)
    if len(fs) > 1:
        for f in fs:
            flow = f.split('/')[-2]
            flow = flow.split('_')[2]
            target = f'{targetdir}/{s}_{flow}_pass.fq.gz'
            f = pathlib.Path(f).absolute().resolve()
            os.symlink(f, target)
    elif t == 'ONT':
        target = f'{targetdir}/{s}_pass.fq.gz'
        [f] = fs
        f = pathlib.Path(f).absolute().resolve()
        os.symlink(f, target)
    elif t == 'ILM':
        for ix,f in enumerate(fs[0]):
            target = f'{targetdir}/{s}.pair.{ix+1}.fq.gz'
            f = pathlib.Path(f).absolute().resolve()
            os.symlink(f, target)
    else:
        raise NotImplementedError("??")

