#!/usr/bin/env bash
# Cluster 100AA smORFs at 90% sequence identity using CD-HIT to generate 90AA smORFs
# Download cdhit: https://github.com/weizhongli/cdhit
#
set -ve
cd-hit \
    -i /work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER/SHD_SMORF_resource/100AA_SmORFs_sequences.faa.gz \
    -o /work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER/clustered_SmORFs \
    -c 0.9 \
    -n 5 \
    -M 16000 \
    -T 16 \
    -d 0
