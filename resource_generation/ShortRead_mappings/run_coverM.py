import pathlib
import tempfile
import subprocess
from jug import utils, TaskGenerator

@TaskGenerator
def create_genomes_directory(include_nces):
    from SemiBin.fasta import fasta_iter
    import pandas as pd
    import gzip
    import os
    odir = ('genomes+nces' if include_nces else 'genomes')
    os.makedirs(odir, exist_ok=True)
    target_dir = pathlib.Path(odir)
    base = pathlib.Path('../../data/ShanghaiDogsMAGs/')
    shd_meta = pd.read_csv('../../data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv', index_col=0)
    for p in base.iterdir():
        if p.suffixes == ['.fna', '.gz']:
            if shd_meta.loc[p.name, 'Representative'] == 'Yes':
                (target_dir / p.name).symlink_to('..' / p)
    if include_nces:
        for h, seq in fasta_iter('../../data/SHD1_NC.fna.gz'):
            with gzip.open(f'{odir}/{h}.fna.gz', 'wt') as out:
                out.write(f'>{h}\n{seq}\n')
    return odir


@TaskGenerator
def run_coverM(b, genome_nce_dir, genome_dir):
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
        r = {}
        for m in ['relative_abundance', 'trimmed_mean', 'covered_fraction']:
            with open(ttable, 'wb') as out:
                ref_dir = (genome_dir if m == 'relative_abundance' else genome_nce_dir)
                subprocess.check_call(
                        ['coverm', 'genome',
                         '-b', sorted_bam,
                         '-d', ref_dir,
                         '-x', 'gz',
                         '-m', m,
                         ],
                        stdout=out)
            t = pd.read_csv(ttable, sep='\t', index_col=0).squeeze()
            t.name = t.name.replace('mapped.sorted','')
            r[m] = t
        return pd.DataFrame(r).fillna(0.)


@TaskGenerator
def collate(cs):
    import pandas as pd
    r = []
    for s,c in cs.items():
        c.columns = pd.MultiIndex.from_product([[s], c.columns])
        r.append(c)
    r = pd.concat(r,axis=1)
    for m in ['relative_abundance', 'trimmed_mean', 'covered_fraction']:
        r.xs(m, axis=1, level=1).to_csv(f'../../intermediate-outputs/external_datasets_mappings/SHD_{m}.tsv.gz', sep='\t')
    return r

bams = utils.cached_glob(f'outputs/mapped_sp_nc/*.bam')

genome_nce_dir = create_genomes_directory(True)
genome_dir = create_genomes_directory(False)
coverages = {}
for b in bams:
    sample = b.split('/')[-1].removesuffix('mapped.bam')
    coverages[sample] = run_coverM(b, genome_nce_dir, genome_dir)

final = collate(coverages)
