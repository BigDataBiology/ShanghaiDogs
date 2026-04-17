from os import makedirs
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt
import deltavis
import dataclasses

plt.rcParams['svg.fonttype'] = 'none'

makedirs('outputs/figures', exist_ok=True)

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
    ax.set_yticklabels(["0", "1M"])
    ax.set_xticks([0, 1_000_000, 2_000_000])
    ax.set_xticklabels(["0", "1M", "2M"])
    ax.set_xlabel(None)
sns.despine(fig)
fig.tight_layout()
fig.savefig('outputs/figures/dotplots8.svg')

fig, ax = plt.subplots(figsize=(7 / IN_2_CM, 8 / IN_2_CM), sharex=True)
ax.clear()
other_genomes = OTHER_GENOMES[0]
delta = deltavis.parse_delta(f'outputs/nucmer/{REF}_vs_{other_genome}.delta')
_, auto_flip = deltavis.auto_orient(delta)
auto_rot = find_rotation(delta)

deltavis.plot_dotplot(delta, ax=ax, flip=auto_flip, rotate=auto_rot)
ax.set_ylabel(f'SHD1_{other_genome}',  fontsize=9)
ax.set_yticks([(+auto_rot) % delta.sections[0].query_len, (1_000_000 + auto_rot ) % delta.sections[0].query_len])
ax.set_yticklabels(["0", "1M"])
ax.set_xticks([0, 1_000_000, 2_000_000])
ax.set_xticklabels(["0", "1M", "2M"])
ax.set_xlabel('SHD1_0457 (bp)', fontsize=9)
sns.despine(fig)
fig.tight_layout()
fig.savefig('outputs/figures/dotplot1.svg')


fig, ax = plt.subplots(figsize=(6 / IN_2_CM, 6 / IN_2_CM), sharex=True)
other_genome = "ECE.113"
ax.clear()
delta = deltavis.parse_delta(f'outputs/nucmer/{REF}_vs_{other_genome}.delta')
_, auto_flip = deltavis.auto_orient(delta)
auto_rot = find_rotation(delta)

deltavis.plot_dotplot(delta, ax=ax, flip=auto_flip, rotate=auto_rot)
ax.set_ylabel(f'SHD1_{other_genome}',  fontsize=9)
ax.set_yticks([(+auto_rot) % delta.sections[0].query_len, (200_000 + auto_rot ) % delta.sections[0].query_len])
ax.set_yticks([0, delta.sections[0].query_len])
ax.set_yticklabels(["0", "0.2M"])
ax.set_xticks([0, 1_000_000, 2_000_000])
ax.set_xticklabels(["0", "1M", "2M"])
ax.set_xlabel(None)
ax.set_xlabel('SHD1_0457 (bp)', fontsize=9)
sns.despine(fig)
fig.tight_layout()
fig.savefig('outputs/figures/dotplot_ec113.svg')

fig, axes = plt.subplots(2,1, figsize=(6 / IN_2_CM, 8 / IN_2_CM), sharex=True)
for ax, other_genome in zip(axes.flat, [OTHER_GENOMES[0], "ECE.113"]):
    ax.clear()
    delta = deltavis.parse_delta(f'outputs/nucmer/{REF}_vs_{other_genome}.delta')
    _, auto_flip = deltavis.auto_orient(delta)
    auto_rot = find_rotation(delta)

    deltavis.plot_dotplot(delta, ax=ax, flip=auto_flip, rotate=auto_rot)
    ax.set_ylabel(f'SHD1_{other_genome}',  fontsize=9)
    if delta.sections[0].query_len > 1_000_000:
        ax.set_yticks([(+auto_rot) % delta.sections[0].query_len, (1_000_000 + auto_rot ) % delta.sections[0].query_len])
        ax.set_yticklabels(["0", "1M"])
    else:
        ax.set_yticks([(+auto_rot) % delta.sections[0].query_len, (200_000 + auto_rot ) % delta.sections[0].query_len])
        ax.set_yticks([0, delta.sections[0].query_len])
        ax.set_yticklabels(["0", "0.2M"])
    ax.set_xticks([0, 1_000_000, 2_000_000])
    ax.set_xticklabels(["0", "1M", "2M"])
    ax.set_xlabel(None)
ax.set_xlabel('SHD1_0457 (bp)', fontsize=9)
sns.despine(fig)
fig.tight_layout()
fig.savefig('outputs/figures/dotplot2.svg')



SAMPLES = [
     'D012', # D012 first
     'D003',
     'D004',
     'D005',
     'D006',
     'D011',
     'D013',
     'D019',
     'D020',
     'D022',
     'D023',
     'D028',
     'D030',
     'D043',
     'D044',
     ]

WINDOW_SIZE = 1000

def bedgraph_to_binned(path):
    regions = []
    max_end = 0
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] != "SHD1_0457_1":
                continue
            start, end, depth = int(parts[1]), int(parts[2]), float(parts[3])
            regions.append((start, end, depth))
            if end > max_end:
                max_end = end
    cov = np.zeros(max_end, dtype=np.float32)
    for start, end, depth in regions:
        cov[start:end] = depth
    n_bins = max_end // WINDOW_SIZE
    binned = cov[: n_bins * WINDOW_SIZE].reshape(n_bins, WINDOW_SIZE).mean(axis=1)
    positions = np.arange(n_bins) * WINDOW_SIZE / 1e6
    return positions, binned



def plot_coverages(axes, make_binary):
    n_ax = 0
    for sample in SAMPLES:
        bedgraph_path = f'outputs/coverage_competitive/{sample}_ShanghaiDogsMAGsSpecies.filtered_SHD1_0457_1.sorted.bedgraph'
        positions, coverage = bedgraph_to_binned(bedgraph_path)
        if coverage.mean() < 2.0:
             continue
        ax = axes.flat[n_ax]
        ax.clear()
        n_ax += 1
        if make_binary:
            coverage = coverage > 0
        color = '#7570b3' if sample != 'D012' else '#e7298a'
        ax.plot(positions, coverage, linewidth=1.0, label=sample, color=color)
        ax.set_xlabel(None)
        ax.set_xlim(0, 1.2)
        ax.set_ylabel(sample, fontsize=7)
        if make_binary:
            ax.set_ylim(-0.1, 1.1)
            ax.set_yticks([0, 1])
            ax.set_yticklabels(["", ""], fontsize=7)
        # add grey rectangles for the regions of interest
        ax.axvspan(0.08, 0.12, color='grey', alpha=0.3)
        ax.axvspan(0.50, 0.52, color='grey', alpha=0.3)
    axes.flat[-1].tick_params(axis="x", labelsize=9)
    axes.flat[-1].set_xlabel('Genome position on SHD1_0457 (Mbp)', fontsize=11)

fig, axes = plt.subplots(8,1, figsize=(10 / IN_2_CM, 8 / IN_2_CM), sharex=True)
plot_coverages(axes, make_binary=True)
sns.despine(fig)
fig.tight_layout()
fig.savefig('outputs/figures/coverages_binary.svg')


fig, axes = plt.subplots(8,1, figsize=(18 / IN_2_CM, (23-8) / IN_2_CM), sharex=True)
plot_coverages(axes, make_binary=False)
for ax in axes.flat:
    ax.set_yscale('log')
    ax.set_xlim(0, 2.3)
sns.despine(fig)
fig.tight_layout()
fig.savefig('outputs/figures/coverages_full.svg')
