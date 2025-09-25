#!/usr/bin/env python3
import itertools
import pandas as pd
import sys
import os
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter
from pathlib import Path
import time

def create_motif_dataframe(num_nucleotides, file_names):
    nucleotides = 'ACGT'
    motifs = [''.join(p) for p in itertools.product(nucleotides, repeat=num_nucleotides)]
    # Clean file names by removing "_motifs.tsv"
    cleaned_file_names = [os.path.basename(f).replace("_motifs.tsv", "") for f in file_names]
    motif_df = pd.DataFrame(0, index=motifs, columns=cleaned_file_names)
    motif_df.index.name = 'motif'  # Name the first column 'motif'
    return motif_df

def process_file(file_path, motifs_set):
    counts = Counter()
    with open(file_path, 'r') as file:
        for line in file:
            motif = line.strip().upper()
            if motif in motifs_set:
                counts[motif] += 1
    return counts

def count_motifs_in_files(directory_path, num_nucleotides, num_cores):
    file_names = [f for f in os.listdir(directory_path) if f.endswith('.tsv')]
    motifs = [''.join(p) for p in itertools.product('ACGT', repeat=num_nucleotides)]
    motifs_set = set(motifs)
    motif_df = create_motif_dataframe(num_nucleotides, file_names)

    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        future_to_file = {
            executor.submit(process_file, os.path.join(directory_path, file_name), motifs_set): file_name
            for file_name in file_names
        }
        for future in tqdm(as_completed(future_to_file), total=len(file_names), desc="Processing Files"):
            file_name = future_to_file[future]
            counts = future.result()
            for motif, count in counts.items():
                cleaned_file_name = os.path.basename(file_name).replace("_motifs.tsv", "")
                motif_df.at[motif, cleaned_file_name] = count

    return motif_df

# Check if correct number of command line arguments are provided
if len(sys.argv) != 5:
    print("Usage: python script.py <directory_of_sample_files> <motif_length> <num_cores> <output_directory>")
    sys.exit(1)

directory_path = sys.argv[1]
motif_length = int(sys.argv[2])
num_cores = int(sys.argv[3])
output_directory = sys.argv[4]

# Create the output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Start timing the script
start_time = time.time()

# Count motifs and create the dataframe
motif_counts_df = count_motifs_in_files(directory_path, motif_length, num_cores)

# End timing the script
end_time = time.time()

# Set the output file names
output_counts_file = os.path.join(output_directory, "motif_counts.tsv")
output_relfreq_file = os.path.join(output_directory, "motifs_rel_freq_wide.tsv")

# Save counts (motifs x libraries)
motif_counts_df.to_csv(output_counts_file, sep='\t')

# ---- NEW: compute and save relative frequencies (libraries x motifs) ----
# Divide each library column by its column sum, then transpose
col_sums = motif_counts_df.sum(axis=0).replace(0, pd.NA)
motifs_rel_freq_wide = motif_counts_df.divide(col_sums, axis=1).T
motifs_rel_freq_wide.index.name = "library"
motifs_rel_freq_wide = motifs_rel_freq_wide.sort_index()
motifs_rel_freq_wide.to_csv(output_relfreq_file, sep='\t')

# Print out the time taken and output locations
print(f"Processed in {end_time - start_time:.2f} seconds.")
print(f"Motif counts saved to {output_counts_file}")
print(f"Relative frequencies (wide) saved to {output_relfreq_file}")
