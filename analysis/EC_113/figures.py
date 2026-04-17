from os import makedirs
import numpy as np
from matplotlib import pyplot as plt
import deltavis
import dataclasses

makedirs('plots', exist_ok=True)

IN_2_CM = 2.54

REF = '0457'
GENOME_WITH_3_CONTIGS = '0522'
OTHER_GENOMES = [
    '0607',
    '0683',
    '0730',
    '0734',
    '0762',
    '0763',
    '0808',
    '0809',
    ]



def rotate_alignments(alignments, rotate, total_query):
    rot_off = rotate % total_query if total_query > 0 else 0
    return [
        dataclasses.replace(aln,
            query_start=(aln.query_start - rot_off) % total_query,
            query_end=(aln.query_end - rot_off) % total_query,
        ) for aln in alignments]

def score_alignments(alignments, rotate, total_len):
    rot_off = rotate % total_len if total_len > 0 else 0
    score = 0.0
    for ix, aln in enumerate(alignments):
        aln_len = np.abs(aln.query_end - aln.query_start)
        aln_qstart = (aln.query_start - rot_off) % total_len
        for other_aln in alignments[:ix]:
            other_len = np.abs(other_aln.query_end - other_aln.query_start)
            other_qstart = (other_aln.query_start - rot_off) % total_len
            monotonic = (1 if aln_qstart > other_qstart else -1)
            score += monotonic * aln_len / total_len * other_len / total_len
    return score

def find_rotation(delta):
    alignments = delta.sections[0].alignments
    lens = np.array(
            [(aln.query_end - aln.query_start) for aln in alignments])
    lens = np.abs(lens)
    total_len = lens.sum()
    lens.sort()
    min_len = lens[np.sum( np.cumsum(lens) < 0.1 * total_len )]
    alignments = [aln for aln in alignments if np.abs(aln.query_end - aln.query_start) >= min_len]
    max_score = 0
    best_rot = 0
    for rot in range(0, total_len, 10_000):
        score = score_alignments(alignments, rot, total_len)
        if np.abs(score) > max_score:
            max_score = np.abs(score)
            best_rot = rot
    return best_rot

fig, axes = plt.subplots(2, 4)
for ax, other_genome in zip(axes.flat, OTHER_GENOMES):
    ax.clear()
    delta = deltavis.parse_delta(f'outputs/nucmer/{REF}_vs_{other_genome}.delta')
    auto_rot, auto_flip = deltavis.auto_orient(delta)
    auto_rot = find_rotation(delta)


    deltavis.plot_dotplot(delta, ax=ax, flip=auto_flip, rotate=auto_rot)
    ax.set_title(other_genome)
    ax.set_ylabel(None)
    #ax.set_yticks([])
    ax.set_xlabel(None)
fig.tight_layout()

