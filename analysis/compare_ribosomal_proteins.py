from jug import TaskGenerator
from jug.utils import cached_glob

RIBOSOMAL_DIR = '../intermediate-outputs/07_ribosomal_genes/barrnap_fasta/'

def compare_seqs(seq1, seq2):
    '''Align two sequences and return average identity'''
    import skbio.alignment
    import skbio.sequence
    aln, score, pos = skbio.alignment.local_pairwise_align_nucleotide(
                seq1=seq1, seq2=seq2, gap_open_penalty=5, gap_extend_penalty=2)
    same = 0
    n = 0
    for pos in aln.iter_positions():
        pos = str(pos)
        if len(pos) != 2:
            print("??")
            break
        if pos[0] == pos[1]:
            same += 1
        n += 1
    return same/n

@TaskGenerator
def compare_all_pairwise(ifile, name_filter):
    '''Compare all 16s sequences in a file pairwise

    Returns a matrix of identities
    '''
    import numpy as np
    from fasta import fasta_iter
    import skbio.sequence
    selected = []
    for name, seq in fasta_iter(ifile):
        if name.startswith(name_filter):
            selected.append(seq)
    selected = [skbio.sequence.DNA(seq) for seq in selected]

    identity = np.zeros((len(selected), len(selected)), dtype=float)
    for i in range(len(selected)):
        for j in range(i+1, len(selected)):
            identity[i,j] = compare_seqs(selected[i], selected[j])
    identity += identity.T

    return identity

for ifile in cached_glob(f'{RIBOSOMAL_DIR}/*/*_ribosomal.fa'):
    compare_all_pairwise(ifile, '16S_rRNA')
    compare_all_pairwise(ifile, '23S_rRNA')
    compare_all_pairwise(ifile,  '5S_rRNA')

