import fasta
from jug import TaskGenerator

@TaskGenerator
def generate_contigs_binned_table():
    import pandas as pd
    metadata = pd.read_csv('../../data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv')

    contig_data = []
    for si in range(0, 53):
        sample = f'D{si:03}'
        if sample == 'D009':
            continue

        sel = metadata.query('Sample == @sample')
        binned = set()
        for b in sel['Bin ID']:
            for h,_ in fasta.fasta_iter(f'../../data/ShanghaiDogsMAGs/{b}', full_header=True):
                _,oh,_,s = h.split(' ')
                assert s == sample
                binned.add(oh)
        for h,seq in fasta.fasta_iter(f'../../data/ShanghaiDogsAssemblies/{sample}_PP1_PolcaCorr.fna.gz'):
            contig_data.append((sample, h, len(seq), h in binned))

    contig_data = pd.DataFrame(contig_data, columns=['Sample', 'Header', 'Length', 'Binned'])
    ofile = '../../intermediate-outputs/Tables/SHD_contigs_binned.csv'
    contig_data.to_csv(ofile, index=False)
    return ofile


@TaskGenerator
def semibin_binned_fraction():
    import pandas as pd
    from os import path
    from macrel import fasta
    from glob import glob

    EX_BINS_BASEDIR = '../../external-data/data/Coelho_2018_bins/'

    samples = sorted([path.basename(f) for f in glob(f'{EX_BINS_BASEDIR}/S3N2Bin/*')])
    contig_data = []

    for sample in samples:
        sm_res = pd.read_csv(f'{EX_BINS_BASEDIR}/SemiBin_benchmark/visualization/Results/Real/CheckM_GUNC/dog/multi_sample/{sample}/SemiBin/result.csv', index_col=0)
        sm_res.eval('med_hi_q = Completeness > 50 and Contamination < 10 and `pass.GUNC`', inplace=True)

        binned = set()
        for f in glob(f'{EX_BINS_BASEDIR}/S3N2Bin/{sample}/*.fa'):
            base = f.split('/')[-1].replace('.fa', '')
            if not sm_res.loc[base, 'med_hi_q']:
                continue


            for h, _ in fasta.fasta_iter(f):
                binned.add(h)

        for h,seq in fasta.fasta_iter(f'{EX_BINS_BASEDIR}/Coelho_2018_dog/{sample}-assembled.fa.gz'):
            contig_data.append((sample, h, len(seq), h in binned))

    contig_data = pd.DataFrame(contig_data, columns=['Sample', 'Header', 'Length', 'Binned'])
    ofile = '../../intermediate-outputs/Tables/Coelho_2018_contigs_semibin_binned.csv'
    contig_data.to_csv(ofile, index=False)
    return ofile

generate_contigs_binned_table()
semibin_binned_fraction()

