import gzip
import os
import csv
import fasta

#paths
base_dir = "../../intermediate-outputs/GMSC_MAPPER"
origins_file = os.path.join(base_dir, "SHD_SMORF_resource", "100AA_SmORFs_origins.tsv.gz")
sequences_file = os.path.join(base_dir, "SHD_SMORF_resource", "100AA_SmORFs_sequences.faa.gz")

def parse_predicted_fasta(file_path):
    metadata_dict = {}
    for header, seq in fasta.fasta_iter(file_path, full_header=True):
        parts = header.split("#")
        if len(parts) >= 5:
            smorf_id = parts[0].strip()
            contig = parts[1].strip()
            start = parts[2].strip()
            end = parts[3].strip()
            strand = "+" if parts[4].strip() == "1" else "-"
            coords = f"{start}-{end}"
            metadata_dict[(smorf_id, contig)] = {
                "contig": contig,
                "coords": coords,
                "strand": strand,
                "sequence": seq.replace(" ", "").replace("\n", "").upper()
            }
    return metadata_dict

def verify_entry(smorf_id, sample_id, contig, coords, strand, sequences, entry_num):
    # Log entry details
    print(f"Entry {entry_num}: SmORF ID: {smorf_id}, Sample ID: {sample_id}, Contig: {contig}, Coords: {coords}, Strand: {strand}")

    if smorf_id not in sequences:
        print(f"Entry {entry_num}: SmORF ID {smorf_id} not found in {sequences_file}")
        return
    sequence = sequences[smorf_id]
    print(f"Entry {entry_num}: Sequence: {sequence[:20]}...")

    # Check sample directory
    sample_dir = os.path.join(base_dir, f"{sample_id}_PP1_PolcaCorr")
    if not os.path.isdir(sample_dir):
        print(f"Entry {entry_num}: Sample directory {sample_dir} does not exist")
        return

    # Find mapped SmORF ID
    mapped_file = os.path.join(sample_dir, "mapped.smorfs.faa")
    if not os.path.exists(mapped_file):
        print(f"Entry {entry_num}: File {mapped_file} does not exist")
        return
    mapped_smorf_id = None
    for header, seq in fasta.fasta_iter(mapped_file):
        if seq.replace(" ", "").replace("\n", "").upper() == sequence:
            mapped_smorf_id = header.split()[0].lstrip(">")
            break
    if not mapped_smorf_id:
        print(f"Entry {entry_num}: No matching sequence for {smorf_id} in {mapped_file}")
        return
    print(f"Entry {entry_num}: Matched SmORF ID: {mapped_smorf_id}")

    # Get predicted file data
    predicted_file = os.path.join(sample_dir, "predicted.filterd.smorf.faa")
    if not os.path.exists(predicted_file):
        print(f"Entry {entry_num}: File {predicted_file} does not exist")
        return
    predicted_metadata = parse_predicted_fasta(predicted_file)
    key = (mapped_smorf_id, contig)
    if key not in predicted_metadata:
        print(f"Entry {entry_num}: SmORF ID {mapped_smorf_id} with contig {contig} not found in {predicted_file}")
        return
    predicted_info = predicted_metadata[key]
    print(f"Entry {entry_num}: Metadata for {mapped_smorf_id} (Contig: {contig}): Coords={predicted_info['coords']}, Strand={predicted_info['strand']}")

    # Verify the data
    if (predicted_info["contig"] == contig and 
        predicted_info["coords"] == coords and 
        predicted_info["strand"] == strand):
        print(f"Entry {entry_num}: Verification successful for {smorf_id} (Mapped SmORF ID: {mapped_smorf_id}):")
        print(f"  Origins: Contig: {contig}, Coords: {coords}, Strand: {strand}")
        print(f"  Predicted: Contig: {predicted_info['contig']}, Coords: {predicted_info['coords']}, Strand: {predicted_info['strand']}")
    else:
        print(f"Entry {entry_num}: Verification failed for {smorf_id} (Mapped SmORF ID: {mapped_smorf_id}):")
        print(f"  Origins: Contig: {contig}, Coords: {coords}, Strand: {strand}")
        print(f"  Predicted: Contig: {predicted_info['contig']}, Coords: {predicted_info['coords']}, Strand: {predicted_info['strand']}")

def verify_last_ten_entries():
    print("Starting QC check for last 10 entries")
    sequences = {}
    for header, seq in fasta.fasta_iter(sequences_file):
        sequences[header.split()[0].lstrip(">")] = seq.replace(" ", "").replace("\n", "").upper()
    if not sequences:
        print(f"No sequences in {sequences_file}")
        return

    # Read origins file
    entries = []
    with gzip.open(origins_file, "rt", newline='') as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader) 
        entries = list(reader)

    if not entries:
        print("Origins file is empty")
        return

    # Process few entries
    last_ten = entries[-10:] if len(entries) >= 10 else entries
    start_entry_num = len(entries) - len(last_ten) + 1

    for idx, entry in enumerate(last_ten, start=start_entry_num):
        if not entry:
            print(f"Entry {idx}: No valid data")
            continue
        smorf_id, sample_id, contig, coords, strand = entry
        verify_entry(smorf_id, sample_id, contig, coords, strand, sequences, idx)
    print("Completed QC check for last 10 entries")

if __name__ == "__main__":
    try:
        verify_last_ten_entries()
    except Exception as e:
        print(f"Error: {str(e)}")
