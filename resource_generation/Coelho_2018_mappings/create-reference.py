import subprocess
import gzip
import os
from glob import glob

basedir = '../../data/ShanghaiDogsMAGs/'
outdir = '../../intermediate-outputs/Coelho_2018_mappings/'
ofile = outdir + 'ShanghaiDogsMAGs.fna.gz'

if os.path.exists(ofile):
    print(f'{ofile} already exists. Exiting.')
    exit()

os.makedirs(outdir, exist_ok=True)
with open(ofile + '.tmp', 'wb') as tmp_out:
    gz_p = subprocess.Popen(['gzip'], stdin=subprocess.PIPE, stdout=tmp_out)
    for file in sorted(glob(basedir + '*.fna.gz')):
        print(f'Processing {file}')
        with gzip.open(file, 'rb') as f:
            while ch := f.read(8192):
                gz_p.stdin.write(ch)
    gz_p.stdin.close()
    if gz_p.wait() != 0:
        raise OSError(f'Error creating {ofile}')

os.rename(ofile + '.tmp', ofile)

