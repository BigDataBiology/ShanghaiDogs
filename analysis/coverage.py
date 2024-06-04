from jug import TaskGenerator
from jug.utils import cached_glob
from collections import defaultdict
from fasta import fasta_iter
from glob import glob
import hashlib
import subprocess

import pandas as pd
import numpy as np
MAP_FA = '/data/anna/animal_metagenome/long-mg-dog/05_dereplication/01_drep/ANI_95/drep_95_MAGs_rf.fa'

def rmean(vs):
    vs.sort()
    n = len(vs)
    return vs[n//10:-n//10].mean()

def hash_string(s):
    sha = hashlib.sha256()
    sha.update(s.encode('ascii'))
    return sha.hexdigest()
    
@TaskGenerator
def get_contig_mapping():
    hash2contig = {}
    for f in glob('/data/anna/animal_metagenome/long-mg-dog/05_dereplication/01_drep/ANI_95/dereplicated_genomes/*.fa.gz'):
        for h, seq in fasta_iter(f):
            sha = hash_string(seq)
            assert sha not in hash2contig
            hash2contig[sha] = (f,h)

    new2old = {}
    for nh,seq in fasta_iter(MAP_FA):
        sha = hash_string(seq)
        h = hash2contig[sha]
        new2old[nh] = h
    assert len(new2old)  == len(hash2contig)
    return new2old


@TaskGenerator
def get_coverage(new2old, ibam):

    coverage = {}
    for nh,seq in fasta_iter(MAP_FA):
        coverage[nh] = np.zeros(len(seq), dtype=np.int32)
        
    p_bedtools = subprocess.Popen([
        'bedtools', 'genomecov',
         '-bga',
         '-ibam', ibam],
         text=True,
         stdout=subprocess.PIPE)
    for line in p_bedtools.stdout:
        contig, start, end, cov = line.strip().split('\t')
        start = int(start)
        end = int(end)
        cov = int(cov)
        coverage[contig][start:end] = cov

    if p_bedtools.wait() != 0:
        raise OSError(f'Error running bedtools')

    coverage_by_bin = defaultdict(list)
    for contig,covs in coverage.items():
        bin,_ = new2old[contig]
        bin = bin.split('/')[-1].split('.')[0]
        coverage_by_bin[bin].append(covs)
        
    for bin,covs in coverage_by_bin.items():
        covs = np.concatenate(covs)
        coverage_by_bin[bin] = covs
    return coverage_by_bin



@TaskGenerator
def summarize_coverage(coverage_by_bin, method):
    import pandas as pd
    if method == 'rmean':
        return pd.Series({b:rmean(covs) for b, covs in coverage_by_bin.items()})
    else:
        raise NotImplementedError(f'Unknown method "{method}"')


@TaskGenerator
def save_coverage_table(coverages, method):
    import pandas as pd
    mimag = pd.read_csv('../data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv', )
    mimag['Bin ID'] = mimag['Bin ID'].str.split('.').str[0]
    bin_rename = mimag[['Original ID', 'Bin ID']].set_index('Original ID').squeeze().to_dict()
    coverages = pd.DataFrame(coverages)
    coverages = coverages.rename(index=bin_rename)
    coverages /= coverages.sum(0)
    coverages.to_csv(f'../intermediate-outputs/repbin_coverage_{method}.tsv', sep='\t')

new2old = get_contig_mapping()
coverages = {}
for ibam in cached_glob('../intermediate-outputs/05_dereplication/map_sp_ANI95/map_SR/*.bam'):
    sample = ibam.split('/')[-1].split('.')[0]
    coverages[sample] = summarize_coverage(get_coverage(new2old, ibam), 'rmean')

save_coverage_table(coverages, 'rmean')

