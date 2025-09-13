def select(infile,inputdir,assembly_dir,outdir):
    from fasta import fasta_iter

    with open(infile,'rt') as f:
        mags = {}
        for line in f:
            linelist = line.strip().split(',')
            if linelist[19].startswith('ALL'):
                mags[linelist[0]] = linelist[18]
            
    for key,value in mags.items():
        mag_name = key.split('.')[0]
        outfile = f'{outdir}{mag_name}.fna'
        final_mag = f'{inputdir}{key}'
        assembly = f'{assembly_dir}{value}/assembly.fasta.gz'
        
        contigs = set()
        for h,seq in fasta_iter(final_mag):
            contig_name = h.split(' ')[1].replace('_polypolish','')
            contigs.add(contig_name)

        with open(outfile,'wt') as out:
            for h,seq in fasta_iter(assembly):
                if h in contigs:
                    out.write(f'>{h}\n{seq}\n')

infile = 'data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv'
inputdir = 'data/ShanghaiDogsMAGs/'
assembly_dir = 'intermediate-outputs/polishing_evaluation/Flye_assemblies/ALL/'
outdir = 'intermediate-outputs/polish_evaluation_mags/flye_bins'
select(infile,inputdir,assembly_dir,outdir)
