import pandas as pd
import os
import gzip
from fasta import fasta_iter

RC_DIR = '../intermediate-outputs/circular_contigs_flye/RC'
MAGS_DIR = '../data/ShanghaiDogsMAGs/'


bin2contig = []
for f in sorted(os.listdir(MAGS_DIR)):
    if not f.endswith('.fna.gz'):
        continue
    mid = f.split('.')[0]
    for line in gzip.open(f'{MAGS_DIR}/{f}', 'rt'):
        if line[0] == '>':
            _, contig, _, sample = line.split()
            bin2contig.append((mid, contig, sample))
    print(mid)

bin2contig = pd.DataFrame(bin2contig, columns=['Bin ID', 'Contig', 'Sample'])
bin2contig['Contig'] = bin2contig['Contig'].str.rsplit('_', n=1).str[0]
contigs_in_mags = set(bin2contig['Sample'] + '_' + bin2contig['Contig'])

high_vs2 = [f.split('.')[0] for f in os.listdir(f'{RC_DIR}/Removed_due_to_full_VS2_score/')]
high_cm2 = [f.split('.')[0] for f in os.listdir(f'{RC_DIR}/Removed_due_to_high_CM2_score/')]

plasmids = pd.read_csv(f'{RC_DIR}/min_1_plasmid_hallmark_gene.tsv', sep='\t', index_col=0)

clusterId2contig = {}
for h, _ in fasta_iter(f'{RC_DIR}/combined_contigs_ref.fasta'):
    cluster_id, contig = h.split('Ω')
    clusterId2contig[cluster_id] = contig

SAMPLE_SPOT_CHECK = 'D006'
contig2seq = {}
data = []
n_cm2 = 0
in_mags = 0
for h, seq in fasta_iter(f'{RC_DIR}/combined_circular.fasta'):
    if h not in clusterId2contig:
        raise ValueError(f'Could not find contig ID for {h}')
    sample, contig = clusterId2contig[h].split('_', 1)

    if h in high_cm2:
        n_cm2 += 1
        continue
    if f'{sample}_{contig}' in contigs_in_mags:
        in_mags += 1
        continue
    if h in plasmids.index:
        cat = 'plasmid'
    elif h in high_vs2:
        cat = 'virus'
    else:
        cat = 'non_categorized'
    data.append([seq, h, sample, contig, cat])
    if sample == SAMPLE_SPOT_CHECK:
        contig2seq[contig] = seq

def sort_key(row):
    if row[4] == 'plasmid':
        cat = 0
    elif row[4] == 'virus':
        cat = 1
    else:
        cat = 2
    return (cat, row[2], row[3])
data.sort(key=sort_key)
with gzip.open('../intermediate-outputs/non_chromosomal/SHD_NC.fna.gz', 'wt') as out:
    for i, row in enumerate(data):
        nh = f'SHD1_NC.{i:03}'
        seq, h, sample, contig, cat = row
        row[0] = nh
        out.write(f'>{nh} {h} {sample} {contig} {cat}\n{seq}\n')
data = pd.DataFrame(data, columns=['Element', 'Working_header', 'Sample', 'Contig', 'Category'])
data.to_csv('../intermediate-outputs/non_chromosomal/non_chromosomal_elements.tsv.gz', sep='\t', index=False)

print(f'Spot check: {SAMPLE_SPOT_CHECK} (contains {len(contig2seq)} elements)')
print('Will check if the sequences are the same and ID matching is correct')
checked = 0
for h,seq in fasta_iter(f'../data/ShanghaiDogsAssemblies/{SAMPLE_SPOT_CHECK}_PP1_PolcaCorr.fna.gz'):
    h, _ = h.rsplit('_', 1)
    if h in contig2seq:
        checked += 1
        assert contig2seq[h] == seq
assert checked == len(contig2seq)
print('Checks passed')
