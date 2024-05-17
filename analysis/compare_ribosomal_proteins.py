from jug import TaskGenerator
from jug import mapreduce
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
def compare_all_pairwise(ifile, name_filter, min_length=None):
    '''Compare all 16s sequences in a file pairwise

    Returns a matrix of identities
    '''
    import numpy as np
    from fasta import fasta_iter
    import skbio.sequence
    selected = []
    for name, seq in fasta_iter(ifile):
        if name.startswith(name_filter):
            if min_length is not None and len(seq) < min_length: continue
            selected.append(seq)
    selected = [skbio.sequence.DNA(seq) for seq in selected]

    identity = np.zeros((len(selected), len(selected)), dtype=float)
    for i in range(len(selected)):
        for j in range(i+1, len(selected)):
            identity[i,j] = compare_seqs(selected[i], selected[j])
    identity += identity.T

    return identity

def min_nz_value(vals):
    from scipy import spatial
    import numpy as np
    if len(vals) <= 1:
        return np.nan
    return spatial.distance.squareform(vals).min()

@TaskGenerator
def build_results_table(fs, mv_16, mv_23, mv_5):
    import pandas as pd
    mags = [f.split('/')[-1] for f in fs]
    return pd.DataFrame({
        'min_id_16s': mv_16,
        'min_id_23s': mv_23,
        'min_id_5s': mv_5,
        }, index=mags)

ribosomal_files = cached_glob(f'{RIBOSOMAL_DIR}/*/*_ribosomal.fa')
ids_16s = []
ids_16s_full = []
ids_23s = []
ids_5s  = []
for ifile in ribosomal_files:
    ids_16s.append(compare_all_pairwise(ifile, '16S_rRNA', min_length=1200))
    ids_16s_full.append(compare_all_pairwise(ifile, '16S_rRNA'))
    ids_23s.append(compare_all_pairwise(ifile, '23S_rRNA'))
    ids_5s .append(compare_all_pairwise(ifile,  '5S_rRNA'))

min_vals_16s = mapreduce.map(min_nz_value, ids_16s, map_step=64)
min_vals_23s = mapreduce.map(min_nz_value, ids_23s, map_step=64)
min_vals_5s  = mapreduce.map(min_nz_value,  ids_5s, map_step=64)

final = build_results_table(ribosomal_files, min_vals_16s, min_vals_23s, min_vals_5s)
