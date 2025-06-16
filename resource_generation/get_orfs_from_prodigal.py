from collections import defaultdict
import gzip
import re
import fasta
from jug import TaskGenerator


BASE_FNA = '../data/ShanghaiDogsAssemblies'
BASE_PRODIGAL = '../intermediate-outputs/Prodigal'
PAT = r'>(\S+)_(\d+) # (\d+) # (\d+) # (-?\d+) # ID=\d+_\d+;partial=(\d+)'

def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


@TaskGenerator
def recover_fna_from_prodigal_faa(sample):
    orfs = defaultdict(list)
    faafile = f'{BASE_PRODIGAL}/{sample}/{sample}_proteins.faa.gz'
    fnafile = f'{BASE_PRODIGAL}/{sample}/{sample}_orfs.fna.gz'
    for line in gzip.open(faafile, 'rt'):
        if line.startswith('>'):
            if m := re.match(PAT, line):
                contig, gene_id, start, end, strand, partial = m.groups()
                orfs[contig].append((f'{contig}_{gene_id}', int(start), int(end), int(strand), partial))
            else:
                raise ValueError(f'Unexpected header format: {line.strip()}')

    contigs = f'{BASE_FNA}/{sample}_PP1_PolcaCorr.fna.gz'
    with gzip.open(fnafile, 'wt') as f_out:
        for h, s in fasta.fasta_iter(contigs):
            for orf_id, start, end, strand, partial in orfs.get(h, []):
                orf_seq = s[start-1:end]
                if strand == -1:
                    orf_seq = rc(orf_seq)
                f_out.write(f'>{orf_id} {start} {end} {strand} {partial}\n{orf_seq}\n')
    return fnafile


@TaskGenerator
def collate_orfs(fnafiles):
    import lib
    ofile = f'{BASE_PRODIGAL}/SHD.ORF.fna.xz'
    ix = 0
    with lib.xz_out(ofile) as out:
        with lib.xz_out(f'{BASE_PRODIGAL}/SHD.ORF.orig.tsv.xz') as o_out:
            o_out.write(f'ORF\tSample\tOriginal_ID\tStart\tEnd\tStrand\tPartial\n'.encode('ascii'))
            for f in fnafiles:
                for h,seq in fasta.fasta_iter(f, full_header=True):
                    orf_id = lib.pad9('SHD.ORF', ix)
                    sample = f.split('/')[-1].split('_')[0]
                    ix += 1
                    orig_h, start, end, strand, partial = h.split()
                    o_out.write(f'{orf_id}\t{sample}\t{orig_h}\t{start}\t{end}\t{strand}\t{partial}\n'.encode('ascii'))
                    out.write(f'>{orf_id}\n{seq}\n'.encode('ascii'))
    return ofile


fnafiles = []
for i in range(0, 53):
    sample = f'D{i:03d}'
    if sample == 'D009': continue
    fnafiles.append(
        recover_fna_from_prodigal_faa(sample))
collate_orfs(fnafiles)
