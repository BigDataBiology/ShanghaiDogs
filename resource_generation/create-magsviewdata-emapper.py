import polars as pl
from collections import Counter
from os import makedirs
BASEDIR = '../intermediate-outputs/09_eggNOG/eggNOG_annot/'
OUTDIR = '../intermediate-outputs/09_eggNOG/eggnog_summary_magsviewdata/'

NR_GENOMES = 2676

makedirs(OUTDIR, exist_ok=True)
for genome_id in range(NR_GENOMES):
    base = f'{BASEDIR}/SHD1_{genome_id:04}'
    gff = pl.read_csv(f'{base}.emapper.genepred.gff',
                    has_header=False,
                    separator='\t',
                    comment_prefix='#',
                    new_columns=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
    gff = gff.rename({'seqid': 'contig'})
    seen = Counter()
    full_geneid = []
    for geneid in gff['contig'].to_list():
        seen[geneid] += 1
        full_geneid.append(f'{geneid}_{seen[geneid]}')
    gff = gff.with_columns(seqid=pl.Series(full_geneid))

    emapper = pl.read_csv(f'{base}.emapper.annotations', has_header=True, separator='\t', comment_prefix='##')
    emapper = emapper.rename({'#query': 'query'})
    assert not set(emapper['query'].to_list()) - set(gff['seqid'].to_list())
    EMAPPER_COLS = ['query', 'COG_category', 'Preferred_name', 'KEGG_ko', 'KEGG_Module']

    data = gff['seqid', 'contig', 'start', 'end', 'strand'].join(
            emapper[EMAPPER_COLS], left_on='seqid', right_on='query')
    data.write_csv(f'{OUTDIR}/SHD1_{genome_id:04}.emapper_summary.tsv', separator='\t')
