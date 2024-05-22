# -*- coding: utf-8 -*-

"""
Created on Tue May 21 16:46:17 2024
@author: Anna Cusco
"""

import os
import re
import shutil

os.chdir('/data/Projects/ShanghaiDogs/')
def process_files(barrnap_file, fasta_file):
    # Read the barrnap data
    barrnap_data = []
    try:
        with open(barrnap_file, 'r') as file:
            for line in file:
                parts = line.strip().split('\t')
                if len(parts) != 10:
                    print(f"Skipping malformed line in barrnap file: {line}")
                    continue
                filename, contig, source, feature, start, end, score, strand, frame, attribute = parts
                attributes = {key: value for key, value in (item.split('=') for item in attribute.split(';'))}
                barrnap_data.append({
                    'contig': contig,
                    'source': source,
                    'feature': feature,
                    'start': int(start),
                    'end': int(end),
                    'strand': strand.strip(),  # Strip whitespace around the strand information
                    'attributes': attributes
                })
        print(f"Successfully read {len(barrnap_data)} entries from barrnap file: {barrnap_file}")
    except Exception as e:
        print(f"Error reading barrnap file: {e}")
        return [], []

    # Read the fasta data
    fasta_data = {}
    try:
        with open(fasta_file, 'r') as file:
            current_key = None
            current_seq = []
            for line in file:
                if line.startswith('>'):
                    if current_key:
                        fasta_data[current_key] = ''.join(current_seq)
                    header = line.strip()
                    match = re.match(r'>(.*)::(.*):(\d+)-(\d+)\((.)\)', header)
                    if match:
                        rRNA_type, contig, start, end, strand = match.groups()
                        if '5S' in rRNA_type:
                            current_key = None  # Skip this entry
                        else:
                            current_key = f"{contig}:{int(start) + 1}-{end}" 
                        current_seq = []
                    else:
                        print(f"Skipping malformed header in fasta file: {header}")
                else:
                    current_seq.append(line.strip())
            if current_key:
                fasta_data[current_key] = ''.join(current_seq)
        print(f"Successfully read {len(fasta_data)} entries from fasta file: {fasta_file}")
    except Exception as e:
        print(f"Error reading fasta file: {e}")
        return [], []

    # Debugging information
    print(f"Read {len(barrnap_data)} entries from barrnap file")
    print(f"Read {len(fasta_data)} entries from fasta file")

    # Map sequences from fasta to barrnap data and categorize them
    partial_sequences = []
    full_sequences = []

    for entry in barrnap_data:
        contig = entry['contig']
        start = entry['start']
        end = entry['end']
        strand = entry['strand']
        attributes = entry['attributes']
        key = f"{contig}:{start}-{end}"
        sequence = fasta_data.get(key)

        if sequence:
            entry['sequence'] = sequence
            header = f">{attributes['Name']}::{contig}:{int(start)-1}-{end}({strand})"
            fasta_entry = f"{header}\n{sequence}"

            if 'partial' in attributes.get('product', '').lower():
                partial_sequences.append(fasta_entry)
            else:
                full_sequences.append(fasta_entry)
        else:
            print(f"No sequence found for key: {key}")

    return partial_sequences, full_sequences

def process_directory(barrnap_directory, fasta_directory):
    for root, dirs, files in os.walk(barrnap_directory):
        for file in files:
            if file.endswith('_barrnap.txt'):
                barrnap_file = os.path.join(root, file)
                print(barrnap_file)
                fasta_file = os.path.join(fasta_directory, file.replace('_barrnap.txt', '_ribosomal.fa'))
                print(fasta_file)

                if os.path.exists(fasta_file):
                    partial_sequences, full_sequences = process_files(barrnap_file, fasta_file)
                    base_filename = os.path.basename(fasta_file).replace('_ribosomal.fa', '')

                    # Save partial sequences
                    if partial_sequences:
                        partial_fasta_file = os.path.join(fasta_directory + '/partial-ribosomal/', base_filename + '_partial.fa')
                        with open(partial_fasta_file, 'w') as partial_file:
                            partial_file.write('\n'.join(partial_sequences))
                        print(f"Created '{partial_fasta_file}' with partial sequences.")

                    # Save full sequences
                    if full_sequences:
                        full_fasta_file = os.path.join(fasta_directory + '/full-ribosomal/', base_filename + '_full.fa')
                        with open(full_fasta_file, 'w') as full_file:
                            full_file.write('\n'.join(full_sequences))
                        print(f"Created '{full_fasta_file}' with full sequences.")
                else:
                    print(f"FASTA file not found for: {barrnap_file}")

# Specify the directories containing the barrnap files and fasta files
barrnap_directory = 'intermediate-outputs/07_ribosomal_genes/barrnap_out/D014'
fasta_directory = 'intermediate-outputs/07_ribosomal_genes/barrnap_fasta/D014'
process_directory(barrnap_directory, fasta_directory)

