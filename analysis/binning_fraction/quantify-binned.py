import fasta
from jug import TaskGenerator
EX_BINS_BASEDIR = '../../external-data/data/Coelho_2018_bins/'

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

@TaskGenerator
def quantify_spire_mags():
    import pandas as pd
    from glob import glob

    spire_meta = pd.read_csv(f'{EX_BINS_BASEDIR}spire_v1_genome_metadata.tsv.gz', index_col=0, sep='\t')

    fas = sorted(glob(f'{EX_BINS_BASEDIR}/Coelho_2018_dog/*-assembled.fa.gz'))
    contig_data = []
    for fa in fas:
        sample = fa.split('/')[-1].split('-')[0]
        binned = set()
        for sp in spire_meta.query('derived_from_sample == @sample').index:
            for h, _ in fasta.fasta_iter(f'{EX_BINS_BASEDIR}/Coelho_2018_dog/{sp}.fa.gz'):
                binned.add(h)
        for h,seq in fasta.fasta_iter(f'{EX_BINS_BASEDIR}/Coelho_2018_dog/{sample}-assembled.fa.gz'):
            contig_data.append((sample, h, len(seq), h in binned))
    contig_data = pd.DataFrame(contig_data, columns=['Sample', 'Header', 'Length', 'Binned'])
    ofile = '../../intermediate-outputs/Tables/Coelho_2018_contigs_spire_binned.csv'
    contig_data.to_csv(ofile, index=False)
    return ofile

@TaskGenerator
def unigene_binned_fraction():
    import pandas as pd
    binned = pd.read_csv('../../intermediate-outputs/Tables/SHD_contigs_binned.csv', )
    binned['SampleHeader'] = binned['Sample'] + '__' + binned['Header']
    binned.set_index('SampleHeader', inplace=True)

    orfs = pd.read_csv('../../intermediate-outputs/Prodigal/SHD.ORF.orig.tsv.xz',
                        sep='\t', names=['ORF', 'Sample', 'Header'])
    orfs['Contig'] = orfs['Header'].str.rsplit('_').str[:-1].str.join('_')
    orfs['SampleHeader'] = orfs['Sample'] + '__' + orfs['Contig']

    orfs = orfs.join(binned[['Binned']], on='SampleHeader')

    clusters = pd.read_csv('../../intermediate-outputs/Prodigal/SHD.clusters.tsv.xz', sep='\t', names=['ORF', 'rel100', 'rep100', 'rel95', 'rep95'])

    orfs.set_index('ORF', inplace=True)
    clusters.set_index('ORF', inplace=True)
    clusters = clusters.join(orfs[['Binned']])

    results = {
            'count' : [
                len(clusters),
                len(set(clusters['rep100'])),
                len(set(clusters['rep95'])),
                ],
            'fraction_binned' : [
                clusters['Binned'].mean(),
                clusters.groupby('rep100')['Binned'].any().mean(),
                clusters.groupby('rep95')['Binned'].any().mean(),
                ],
            }
    return pd.DataFrame(results, index=['ORF', '100NT', '95NT'])

generate_contigs_binned_table()
semibin_binned_fraction()
quantify_spire_mags()
unigene_binned_fraction()
