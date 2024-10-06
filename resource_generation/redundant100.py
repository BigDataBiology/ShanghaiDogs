from macrel import fasta
from collections import defaultdict
import gzip

from lib import pad9, xz_out

seqs = list(fasta.fasta_iter('../intermediate-outputs/Prodigal/SHD.ORF.fna.xz'))
def k(h_s):
    h,s = h_s
    return (len(s), s, h)

print(f'Loaded {len(seqs)}')

seqs.sort(key=k)

deduped = []
matches = []
prev = ''
h_prev = None
for h_seq in seqs:
    h,seq = h_seq
    if seq == prev:
        matches.append((h, '=', h_prev))
    else:
        deduped.append(h_seq)
        prev = seq
        h_prev = h

del seqs

min_len = len(deduped[0][1])

def rolling_hashes(seq, window_size):
    for i in range(len(seq) - window_size + 1):
        yield hash(seq[i:i+window_size])

hash_matches = defaultdict(list)
eliminated = set()

print(f'Deduped: {len(deduped)}')
for ix, (h, seq) in enumerate(deduped):
    if ix % 100_000 == 0:
        print(pad9('Dedupe iteration', ix))

    candidates = set()
    for ha in set(rolling_hashes(seq, min_len)):
       candidates.update(hash_matches[ha])
       hash_matches[ha].append(ix)
    for ix2 in candidates:
        if ix2 in eliminated:
            continue
        h2, seq2 = deduped[ix2]
        assert ix2 < ix
        if seq2 in seq:
            matches.append((h2, 'C', h))
            assert h2 != h
            eliminated.add(ix2)


with xz_out('../intermediate-outputs/Prodigal/SHD.100NT.fna.xz') as f:
    n_ix = 0
    for ix, (h, seq) in enumerate(deduped):
        if ix in eliminated:
            continue
        n_h = pad9('SHD.100NT', n_ix)
        f.write(f'>{n_h} {h}\n{seq}\n'.encode('utf-8'))
        n_ix += 1
with xz_out('../intermediate-outputs/Prodigal/SHD.100NT.matches.xz') as f:
    for h1, c, h2 in matches:
        f.write(f'{h1}\t{c}\t{h2}\n'.encode('utf-8'))
