from os import makedirs
import seaborn as sns
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

def find_rotation(delta):
    alignments = delta.sections[0].alignments
    alignments = np.array([(aln.query_start, aln.query_end - aln.query_start) for aln in alignments])
    lens = np.abs(alignments[:, 1]).copy()
    total_len = lens.sum()
    lens.sort()
    min_len = lens[np.sum( np.cumsum(lens) < 0.1 * total_len )]
    alignments = alignments.astype(np.float64)
    alignments = alignments[np.abs(alignments.T[1]) >= min_len]


    max_score = 0
    best_rot = 0
    alignments.T[1] /= total_len
    for rot in range(0, total_len, 1_000):
        rot_off = rot % total_len if total_len > 0 else 0
        ralignments = alignments.copy()
        ralignments[:, 0] = (ralignments[:, 0] - rot_off) % total_len
        score = 0.0
        for ix in range(len(ralignments)):
            aln = ralignments[ix]
            other_aln = ralignments[:ix]
            monotonic = np.choose(aln[0] > other_aln[:, 0], [1, -1])
            score += np.sum(monotonic * (aln[1] * other_aln[:, 1]))
        if np.abs(score) > max_score:
            max_score = np.abs(score)
            best_rot = rot
    return best_rot

fig, axes = plt.subplots(2, 4, figsize=(18 / IN_2_CM, 8 / IN_2_CM), sharex=True)
assert len(OTHER_GENOMES) == len(axes.flat)
for ax, other_genome in zip(axes.flat, OTHER_GENOMES):
    ax.clear()
    delta = deltavis.parse_delta(f'outputs/nucmer/{REF}_vs_{other_genome}.delta')
    _, auto_flip = deltavis.auto_orient(delta)
    auto_rot = find_rotation(delta)

    deltavis.plot_dotplot(delta, ax=ax, flip=auto_flip, rotate=auto_rot)
    ax.set_ylabel(None)
    ax.set_title(f'SHD1_{other_genome}',  fontsize=9)
    ax.set_yticks([(+auto_rot) % delta.sections[0].query_len, (1_000_000 + auto_rot ) % delta.sections[0].query_len])
    ax.set_yticklabels(["0", "1m"])
    ax.set_xticks([0, 1_000_000, 2_000_000])
    ax.set_xticklabels(["0", "1m", "2m"])
    ax.set_xlabel(None)
sns.despine(fig)
fig.tight_layout()
fig.savefig('plots/dotplots8.svg')

