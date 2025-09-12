import os
from atomicwrites import atomic_write
import pandas as pd
import numpy as np

data_list = ['D000', 'D001', 'D002', 'D003', 'D004', 'D005', 'D006', 'D007', 'D008', 'D010', 'D011', 'D012', 'D013', 'D014',
              'D015', 'D016', 'D017', 'D018', 'D019', 'D020', 'D021', 'D022', 'D023', 'D024', 'D025', 'D026', 'D027', 'D028', 'D029',
              'D030', 'D031', 'D032', 'D033', 'D034', 'D035', 'D036', 'D037', 'D038', 'D039', 'D040', 'D041', 'D042', 'D043', 'D044',
              'D045', 'D046', 'D047', 'D048', 'D049', 'D050', 'D051', 'D052']

def fasta_iter(fname, full_header=False):
    '''Iterate over a (possibly gzipped) FASTA file

    Parameters
    ----------
    fname : str
        Filename.
            If it ends with .gz, gzip format is assumed
            If .bz2 then bzip2 format is assumed
            if .xz, then lzma format is assumerd
    full_header : boolean (optional)
        If True, yields the full header. Otherwise (the default), only the
        first word

    Yields
    ------
    (h,seq): tuple of (str, str)
    '''
    header = None
    chunks = []
    if fname.endswith('.gz'):
        import gzip
        op = gzip.open
    elif fname.endswith('.bz2'):
        import bz2
        op = bz2.open
    elif fname.endswith('.xz'):
        import lzma
        op = lzma.open
    else:
        op = open
    with op(fname, 'rt') as f:
        for line in f:
            if line[0] == '>':
                if header is not None:
                    yield header,''.join(chunks)
                line = line[1:].strip()
                if not line:
                    header = ''
                elif full_header:
                    header = line.strip()
                else:
                    header = line.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, ''.join(chunks)

def generate_split():
    for ref in data_list:
        os.makedirs(f"/data/sjpan/Soil/{ref}", exist_ok=True)
        ref_fasta = f'/data/anna/animal_metagenome/long-mg-dog/02_polishing/02_POLCA/{ref}/{ref}_PP1_PolcaCorr.fasta'
        output_dir = f"/data/sjpan/Soil/{ref}"
        os.makedirs(output_dir, exist_ok=True)
        with atomic_write(f"{output_dir}/contig_split.fa") as ofile:
            for h, seq in fasta_iter(ref_fasta):
                half = len(seq) // 2
                h1 = h + "_1"
                seq1 = seq[:half]
                h2 = h + "_2"
                seq2 = seq[half:]
                ofile.write(f">{h1}\n{seq1}\n")
                ofile.write(f">{h2}\n{seq2}\n")

def mapping():
    for ref in data_list:
        ref_fasta = f'/data/sjpan/Soil/{ref}/contig_split.fa'
        # read1 = f'/data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ngs/{data_list[0]}_350_trim_filter.pair.1.fq.gz'
        # read2 = f'/data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ngs/{data_list[0]}_350_trim_filter.pair.2.fq.gz'
        # os.system(f"/data/sjpan/AEMB/build/strobealign --create-index {ref_fasta} {read1} {read2}")
        os.makedirs(f"/data/sjpan/Soil/{ref}", exist_ok=True)
        for query in data_list:
            read1 = f'/data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ngs/{query}_350_trim_filter.pair.1.fq.gz'
            read2 = f'/data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ngs/{query}_350_trim_filter.pair.2.fq.gz'
            os.system(f"/data/sjpan/AEMB/build/strobealign {ref_fasta} {read1} {read2} -t 64 -R 6 -x > "
                      f"/data/sjpan/Soil/{ref}/{query}.txt")

def generate_csv():
    for ref in data_list:
        print(ref)
        sample_list = []

        for query in data_list:
            sample1 = pd.read_csv(f"/data/sjpan/Soil/{ref}/{query}.txt", sep='\t', index_col=0, header=None)
            sample1.index.name = None
            sample1.columns = [f'sample{query}']
            sample1 = sample1.apply(lambda x: x + 1e-5)
            abun_scale = (sample1.mean() / 100).apply(np.ceil) * 100
            sample1 = sample1.div(abun_scale)
            sample_list.append(sample1)

        sample = pd.concat(sample_list, axis=1)
        sample.to_csv(f'/data/sjpan/Soil/{ref}/cov_split.csv')

        index_name = sample.index.tolist()
        data_index_name = []
        for i in range(len(index_name) // 2):
            data_index_name.append(index_name[2*i][:-2])

        data_split = sample.values

        data = (data_split[::2] + data_split[1::2]) / 2
        columns = [index for index in data_list]
        data = pd.DataFrame(data, index=data_index_name, columns=columns)
        data.to_csv(f'/data/sjpan/Soil/{ref}/cov.csv')

def generate_data():
    for ref in data_list:
        print(ref)
        data = pd.read_csv(f"/data/anna/animal_metagenome/long-mg-dog/04_binning/00_SemiBin2/{ref}/data.csv", index_col = 0)
        data1 = data.iloc[:, :136]

        cov = pd.read_csv(f"/data/sjpan/Soil/{ref}/cov.csv", index_col=0)
        data1 = pd.merge(data1, cov, how='inner', on=None,
                        left_index=True, right_index=True, sort=False, copy=True)
        data1.to_csv(f"/data/sjpan/Soil/{ref}/data.csv")
        print(data.shape, data1.shape)
        assert data.shape[0] == data1.shape[0]
        assert data.shape[1] == data1.shape[1] - 50

        data = pd.read_csv(f'/data/anna/animal_metagenome/long-mg-dog/04_binning/00_SemiBin2/{ref}/data_split.csv', index_col = 0)
        data1 = data.iloc[:, :136]

        cov = pd.read_csv(f'/data/sjpan/Soil/{ref}/cov_split.csv', index_col=0)
        data1 = pd.merge(data1, cov, how='inner', on=None,
                        left_index=True, right_index=True, sort=False, copy=True)
        data1.to_csv(f'/data/sjpan/Soil/{ref}/data_split.csv')
        print(data.shape, data1.shape)
        assert data.shape[0] == data1.shape[0]
        assert data.shape[1] == data1.shape[1] - 52

def generate_training_data():
    os.makedirs(f"/data/sjpan/Soil/training", exist_ok=True)
    for ref in data_list:
        print(ref)
        os.makedirs(f"/data/sjpan/Soil/training/{ref}", exist_ok=True)
        os.system(f"cp /data/sjpan/Soil/{ref}/data.csv /data/sjpan/Soil/training/{ref}")
        os.system(f"cp /data/sjpan/Soil/{ref}/data_split.csv /data/sjpan/Soil/training/{ref}")
        os.system(f"cp /data/anna/animal_metagenome/long-mg-dog/02_polishing/02_POLCA/{ref}/{ref}_PP1_PolcaCorr.fasta /data/sjpan/Soil/training/{ref}")

def run_SemiBin2():
    for ref in data_list:
        print(ref)
        # os.system(f"rm -rf /public/home/pansj/abundance_estimation/training/{ref}/output_recluster_bins")
        # os.system(f"rm -rf /public/home/pansj/abundance_estimation/training/{ref}/output_prerecluster_bins")
        # if os.path.exists(f"/public/home/pansj/abundance_estimation/training/{ref}/checkm2_output"):
        #     os.system(f"rm -rf /public/home/pansj/abundance_estimation/training/{ref}/checkm2_output")
        # os.system(f"SemiBin2 train_self --data "
        #           f"/public/home/pansj/abundance_estimation/training/{ref}/data.csv "
        #           f"--data-split /public/home/pansj/abundance_estimation/training/{ref}/data_split.csv "
        #           f"-o /public/home/pansj/abundance_estimation/training/{ref} -t 4")
        os.system(f"SemiBin2 bin_long "
                  f"-i /public/home/pansj/abundance_estimation/training/{ref}/{ref}_PP1_PolcaCorr.fasta "
                  f"-o /public/home/pansj/abundance_estimation/training/{ref} --data "
                  f"/public/home/pansj/abundance_estimation/training/{ref}/data.csv --model "
                  f"/public/home/pansj/abundance_estimation/training/{ref}/model.h5 "
                  f"--orf-finder fraggenescan -t 4 --write-pre-reclustering-bins --compression none")


def run_checkm2():
    for ref in data_list:
        print(ref)
        os.chdir(f"/public/home/pansj/abundance_estimation/training/{ref}/")
        os.system(f"/public/home/pansj/software/CheckM2-main/bin/checkm2 predict --threads 30 --input output_bins --output-directory checkm2_output -x .fa")

data_list1 = ['D000', 'D001', 'D002', 'D003', 'D004', 'D005', 'D006', 'D007', 'D008', 'D010', 'D011', 'D012', 'D013', 'D014',
              'D015', 'D016', 'D017', 'D018', 'D019', 'D020', 'D021', 'D022', 'D023', 'D024', 'D025', 'D026', 'D027', 'D028', 'D029']

def result_single_sample():
    num_high = 0
    num_medium = 0
    for ref in data_list:
        print(ref)
        SCG = pd.read_csv(f"/data/anna/animal_metagenome/long-mg-dog/04_binning/00_SemiBin2/{ref}/checkm2/quality_report.tsv", index_col=0, sep='\t')
        index_name = SCG.index
        index_name = [str(temp) for temp in index_name]
        SCG.index = index_name
        merge = SCG
        high_qulaity = merge[(merge['Completeness'].astype(float) > float(90)) & (
                merge['Contamination'].astype(float) < float(0.05 * 100))].shape[0]
        medium_quality = merge[(merge['Completeness'].astype(float) >= float(50)) & (merge['Contamination'].astype(float) < float(0.1 * 100))].shape[0] - high_qulaity
        print(high_qulaity, medium_quality)
        num_high += high_qulaity
        num_medium += medium_quality
    print(num_high, num_medium)

def result_multi_sample():
    num_high = 0
    num_medium = 0
    for ref in data_list:
        print(ref)
        SCG = pd.read_csv(f"/public/home/pansj/abundance_estimation/training/{ref}/checkm2_output/quality_report.tsv", index_col=0, sep='\t')
        index_name = SCG.index
        index_name = [str(temp) for temp in index_name]
        SCG.index = index_name
        merge = SCG
        high_qulaity = merge[(merge['Completeness'].astype(float) > float(90)) & (
                merge['Contamination'].astype(float) < float(0.05 * 100))].shape[0]
        medium_quality = merge[(merge['Completeness'].astype(float) >= float(50)) & (merge['Contamination'].astype(float) < float(0.1 * 100))].shape[0] - high_qulaity
        print(high_qulaity, medium_quality)
        num_high += high_qulaity
        num_medium += medium_quality
    print(num_high, num_medium)

if __name__ == '__main__':
    # mapping()
    # generate_split()

    # generate_csv()
    # generate_data()
    # generate_training_data()
    # run_SemiBin2()
    # run_checkm2()
    result_single_sample()
    # result_multi_sample()