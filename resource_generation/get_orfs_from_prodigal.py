from collections import defaultdict
import gzip
import re
import fasta
from jug import TaskGenerator


BASE_FNA = '../data/ShanghaiDogsAssemblies'
BASE_PRODIGAL = '../intermediate-outputs/Prodigal'
PAT = r'>(\S+)_(\d+) # (\d+) # (\d+) # (-?\d+)'
def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

def pad9(ix):
    ix = f'{ix:09}'
    return f'{ix[:3]}_{ix[3:6]}_{ix[6:]}'

@TaskGenerator
def recover_fna_from_prodigal_faa(sample):
    orfs = defaultdict(list)
    faafile = f'{BASE_PRODIGAL}/{sample}/{sample}_proteins.faa.gz'
    fnafile = f'{BASE_PRODIGAL}/{sample}/{sample}_orfs.fna.gz'
    for line in gzip.open(faafile, 'rt'):
        if line.startswith('>'):
            if m := re.match(PAT, line):
                contig, gene_id, start, end, strand = m.groups()
                orfs[contig].append((f'{contig}_{gene_id}', int(start), int(end), int(strand)))

    contigs = f'{BASE_FNA}/{sample}_PP1_PolcaCorr.fna.gz'
    with gzip.open(fnafile, 'wt') as f_out:
        for h, s in fasta.fasta_iter(contigs):
            for orf_id, start, end, strand in orfs.get(h, []):
                orf_seq = s[start-1:end]
                if strand == -1:
                    orf_seq = rc(orf_seq)
                f_out.write(f'>{orf_id}\n{orf_seq}\n')
    return fnafile

@TaskGenerator
def collate_orfs(fnafiles):
    import pandas as pd
    ofile = f'{BASE_PRODIGAL}/SHD_orfs.fna.gz'
    ix = 0
    metadata = []
    with gzip.open(ofile, 'wt') as out:
        for f in fnafiles:
            for h,seq in fasta.fasta_iter(f):
                orf_id = 'SHD_ORF.'+pad9(ix)
                sample = f.split('/')[-1].split('_')[0]
                ix += 1
                metadata.append((orf_id, sample, h))
                out.write(f'>{orf_id}\n{seq}\n')
    metadata = pd.DataFrame(metadata, columns=['orf_id', 'sample', 'orig_id'])
    metadata.to_csv(f'{BASE_PRODIGAL}/SHD_orfs.orig.tsv.gz', sep='\t')
    return ofile


fnafiles = []
for i in range(0, 53):
    sample = f'D{i:03d}'
    if sample == 'D009': continue
    fnafiles.append(
        recover_fna_from_prodigal_faa(sample))
collate_orfs(fnafiles)
