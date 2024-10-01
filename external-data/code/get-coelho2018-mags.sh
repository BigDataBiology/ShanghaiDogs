#!/usr/bin/env bash
set -ev

cd ../data
mkdir -p Coelho_2018_bins
cd Coelho_2018_bins

wget 'https://zenodo.org/records/5181385/files/SemiBin_multi.tar.gz?download=1'
tar -xvzf 'SemiBin_multi.tar.gz?download=1'
rm 'SemiBin_multi.tar.gz?download=1'



wget 'https://swifter.embl.de/~fullam/spire/compiled/Coelho_2018_dog_spire_v1_assemblies.tar'
wget 'https://swifter.embl.de/~fullam/spire/compiled/Coelho_2018_dog_spire_v1_MAGs.tar'

tar -xvf Coelho_2018_dog_spire_v1_assemblies.tar
tar -xvf Coelho_2018_dog_spire_v1_MAGs.tar

rm Coelho_2018_dog_spire_v1_assemblies.tar
rm Coelho_2018_dog_spire_v1_MAGs.tar

wget 'wget https://swifter.embl.de/~fullam/spire/metadata/spire_v1_genome_metadata.tsv.gz'

git clone https://github.com/BigDataBiology/SemiBin_benchmark
cd SemiBin_benchmark/visualization
tar -xvzf Results.tar.gz
