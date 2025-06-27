
### Generate 90AA smORFs based on the 100AA smORFs uses CD - HIT
./cd-hit -i /work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER/Mapped_SmORFs_sequences.faa -o /work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER/clustered_SmORFs -c 0.9 -n 5 -M 16000 -T 16 -d 0
