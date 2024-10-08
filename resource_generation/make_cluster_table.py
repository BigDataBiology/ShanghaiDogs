import pandas as pd
import fasta
import lib

nr_orfs = sum(1 for _ in fasta.fasta_iter('../intermediate-outputs/Prodigal/SHD.ORF.fna.xz'))

m100_to_95 = pd.read_csv('../intermediate-outputs/Prodigal/SHD.95NT.matches.tsv.xz', sep='\t' , header=None, index_col=0)
morf_red = pd.read_csv('../intermediate-outputs/Prodigal/SHD.100NT.matches.xz', sep='\t', header=None, )
morf_to_m100 = []
for h,_ in fasta.fasta_iter('../intermediate-outputs/Prodigal/SHD.100NT.fna.xz', full_header=True):
    h_100,h_orf = h.split()
    morf_to_m100.append((h_orf,h_100))

assert nr_orfs == len(morf_red) + len(morf_to_m100)

assert len(morf_to_m100) == len(set(morf_to_m100))
assert len(morf_to_m100) == len(set(m for m,_ in morf_to_m100))
assert len(morf_to_m100) == len(set(m for _,m in morf_to_m100))
morf_to_m100_d = dict(morf_to_m100)
assert morf_red[0].map(morf_to_m100_d.get).isna().all()

morf_red_d = morf_red.set_index(0).to_dict(orient='index')


with lib.xz_out('../intermediate-outputs/Prodigal/SHD.clusters.tsv.xz') as out:
    for ix in range(nr_orfs):
        ix = lib.pad9('SHD.ORF', ix)
        cur = ix
        rel100 = '='
        while cur not in morf_to_m100_d:
            step = morf_red_d[cur]
            if step[1] == 'C': rel100 = 'C'
            cur = step[2]
        rep100 = morf_to_m100_d[cur]
        rel95,rep95 = m100_to_95.loc[rep100]
        out.write(f'{ix}\t{rel100}\t{rep100}\t{rel95}\t{rep95}\n'.encode('ascii'))
