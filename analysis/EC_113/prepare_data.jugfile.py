from jug import TaskGenerator, Task
GENOMES = [
    '0457',
    '0522',
    '0607',
    '0683',
    '0730',
    '0734',
    '0762',
    '0763',
    '0808',
    '0809',
    '1362',
    '1670',
    '2109',
    '2260',
    '2290',
    ]

@TaskGenerator
def retrieve_genome(gid):
    import os
    import gzip
    os.makedirs('data/', exist_ok=True)
    dest = f'data/SHD1_{gid}.fna'
    original = f'../../data/ShanghaiDogsMAGs/SHD1_{gid}.fna.gz'
    with open(dest, 'wb') as f:
        with gzip.open(original, 'rb') as g:
            while chunk := g.read(8192):
                f.write(chunk)
    return dest

@Task
def retrieve_ec_113():
    from fasta import fasta_iter
    dest = 'data/SHD1_EC.113.fna'
    with open('data/SHD1_EC.113.fna', 'w') as f:
        for h, seq in fasta_iter('../../data/ShanghaiDogs_OtherResources/SHD1_EC.fna.gz'):
            if h == 'SHD1_EC.113':
                f.write(f'>{h}\n{seq}\n')
    return dest


for g in GENOMES:
    retrieve_genome(g)

