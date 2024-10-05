import subprocess
import lzma
import fasta

import tempfile
from lib import pad9, xz_out


with tempfile.TemporaryDirectory() as tdir:
    fnaf = tdir + '/SHD.100NT.fna'
    with lzma.open('../intermediate-outputs/Prodigal/SHD.100NT.fna.xz', 'rb') as ifile, \
            open(fnaf, 'wb') as ofile:
                while ch := ifile.read(8192):
                    ofile.write(ch)
    subprocess.check_call(
        ['cd-hit-est',
            '-c', '0.95',
            '-d', '0',
            '-g', '1',
            '-G', '0',
            '-aS', '0.9',
            '-T', '64',
            '-i', fnaf,
            '-M', '256000',
            '-o', '../intermediate-outputs/Prodigal/cd_hit_out'
         ])

with xz_out('../intermediate-outputs/Prodigal/SHD.95NT.fna.xz') as out:
    for n_ix, (h,seq) in enumerate(fasta.fasta_iter('../intermediate-outputs/Prodigal/cd_hit_out')):
        n_h = pad9('SHD.95NT', n_ix)
        out.write(f'>{n_h}\n{seq}\n'.encode('ascii'))

def iter_clusters(ifile):
    cur = None
    with open(ifile) as f:
        for line in f:
            if line[0] == '>':
                if cur:
                    yield cur
                cur = []
            else:
                _,_,gid,*rest = line.split()
                gid = gid[1:-3]
                cur.append((gid, rest[0] == '*'))
    if cur:
        yield cur

with xz_out('../intermediate-outputs/Prodigal/SHD.95NT.matches.tsv.xz') as out:
    for ix, cl in enumerate(iter_clusters('../intermediate-outputs/Prodigal/cd_hit_out.clstr')):
        [rep] = [g for g,is_r in cl if is_r]
        rep95 = pad9('SHD.95NT', n_ix - ix)
        out.write(f'{rep}\t=\t{rep95}\n'.encode('ascii'))
        for g,is_r in cl:
            if is_r: continue
            out.write(f'{g}\tR\t{rep95}\n'.encode('ascii'))

