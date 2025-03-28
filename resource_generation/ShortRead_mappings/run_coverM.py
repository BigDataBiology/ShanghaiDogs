import pathlib
import tempfile
import subprocess
from jug import utils, TaskGenerator

@TaskGenerator
def create_genomes_directory():
    from SemiBin.fasta import fasta_iter
    import pandas as pd
    import gzip
    import os
    odir = 'genomes'
    os.makedirs(odir, exist_ok=True)
    target_dir = pathlib.Path(odir)
    base = pathlib.Path('../../data/ShanghaiDogsMAGs/')
    shd_meta = pd.read_csv('../../data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv', index_col=0)
    for p in base.iterdir():
        if p.suffixes == ['.fna', '.gz']:
            if shd_meta.loc[p.name, 'Representative'] == 'Yes':
                (target_dir / p.name).symlink_to('..' / p)
    for h, seq in fasta_iter('../../data/SHD1_NC.fna.gz'):
        with gzip.open(f'{odir}/{h}.fna.gz', 'wt') as out:
            out.write(f'>{h}\n{seq}\n')
    return odir


@TaskGenerator
def run_coverM(b, genome_dir):
    import pandas as pd
    with tempfile.TemporaryDirectory() as tdir:
        b = pathlib.PurePath(b)
        sorted_bam = f'{tdir}/{b.stem}.sorted.bam'
        samtools_c = ['samtools', 'sort', b]
        with open(sorted_bam, 'wb') as sbam:
            subprocess.check_call(
                    ['samtools', 'sort', b],
                    stdout=sbam)

        ttable = f'{tdir}/cover.tsv'
        with open(ttable, 'wb') as out:
            subprocess.check_call(
                    ['coverm', 'genome',
                     '-b', sorted_bam,
                     '-d', genome_dir,
                     '-x', 'gz'],
                    stdout=out)
        r = pd.read_csv(ttable, sep='\t', index_col=0).squeeze()
        r.name = r.name.replace('mapped.sorted','')
        return r


@TaskGenerator
def collate(cs):
    import pandas as pd
    pd.DataFrame(cs).to_csv('../../data/ShanghaiDogsTables/SHD_relative_abundances.tsv.gz', sep='\t')

bams = utils.cached_glob(f'outputs/mapped_sp_nc/*.bam')

coverages = []
for b in bams:
    coverages.append(run_coverM(b, create_genomes_directory()))

final = collate(coverages)
