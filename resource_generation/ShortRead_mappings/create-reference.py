import pandas as pd
import subprocess
import gzip
import os
from glob import glob
from jug import TaskGenerator

basedir = '../../data/ShanghaiDogsMAGs/'
outdir = '../../intermediate-outputs/Coelho_2018_mappings/'
plasmid_data = '../../data/SHD1_NC.fna.gz'
shd_meta = pd.read_csv('../../data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv', index_col=0)

@TaskGenerator
def create_reference(ofile, only_reps=False, include_nonchrom=False):
    if os.path.exists(ofile):
        print(f'{ofile} already exists. Exiting.')
        return ofile

    os.makedirs(outdir, exist_ok=True)

    with open(ofile + '.tmp', 'wb') as tmp_out:
        gz_p = subprocess.Popen(['gzip'], stdin=subprocess.PIPE, stdout=tmp_out)
        for file in sorted(glob(basedir + '*.fna.gz')):
            print(f'Processing {file}')
            bin_name = os.path.basename(file)
            if only_reps and shd_meta.loc[bin_name, 'Representative'] != 'Yes':
                print(f'Skipping {bin_name} as it is not a representative.')
                continue
            with gzip.open(file, 'rb') as f:
                while ch := f.read(8192):
                    gz_p.stdin.write(ch)
        if include_nonchrom:
            print("Including plasmids...")
            with gzip.open(plasmid_data, 'rb') as f:
                while ch := f.read(8192):
                    gz_p.stdin.write(ch)
        gz_p.stdin.close()
        if gz_p.wait() != 0:
            raise OSError(f'Error creating {ofile}')

    os.rename(ofile + '.tmp', ofile)
    return ofile

create_reference(outdir + 'ShanghaiDogsMAGs.fna.gz', only_reps=False)
create_reference(outdir + 'ShanghaiDogsMAGs+NonChrom.fna.gz', only_reps=False, include_nonchrom=True)
create_reference(outdir + 'ShanghaiDogsMAGsSpecies+NonChrom.fna.gz', only_reps=True, include_nonchrom=True)
create_reference(outdir + 'ShanghaiDogsMAGsSpecies.fna.gz', only_reps=True)

