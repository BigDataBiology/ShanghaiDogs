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

generate_contigs_binned_table()
