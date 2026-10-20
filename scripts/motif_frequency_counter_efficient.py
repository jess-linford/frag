#!/usr/bin/env python3

"""Combine per-sample 5' motif count tables into count and frequency matrices.

Input files are the output of extract_5p_motifs.sh and have two tab-separated
columns with no header:

    motif    count

The special OTHER row is ignored when constructing the A/C/G/T motif matrix.
"""

import itertools
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from tqdm import tqdm


def all_motifs(motif_length):
    """Return all canonical A/C/G/T motifs of the requested length."""
    return [
        "".join(bases)
        for bases in itertools.product("ACGT", repeat=motif_length)
    ]


def sample_name(file_path):
    """Convert SAMPLE_motifs.tsv to SAMPLE."""
    return Path(file_path).name.removesuffix("_motifs.tsv")


def read_motif_counts(file_path, motifs_set):
    """Read one motif<TAB>count file and return canonical motif counts."""
    counts = {}

    with open(file_path, "r") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) != 2:
                raise ValueError(
                    f"{file_path}:{line_number}: expected 2 tab-separated "
                    f"fields (motif and count), found {len(fields)}"
                )

            motif, count_text = fields
            motif = motif.upper()

            # extract_5p_motifs.sh reports non-canonical motifs as OTHER.
            # Those reads are intentionally excluded from the canonical motif
            # count and relative-frequency matrices.
            if motif == "OTHER":
                continue

            if motif not in motifs_set:
                raise ValueError(
                    f"{file_path}:{line_number}: unexpected motif {motif!r}"
                )

            try:
                count = int(count_text)
            except ValueError as exc:
                raise ValueError(
                    f"{file_path}:{line_number}: invalid count {count_text!r}"
                ) from exc

            if count < 0:
                raise ValueError(
                    f"{file_path}:{line_number}: count cannot be negative"
                )

            counts[motif] = count

    return counts


def combine_motif_counts(directory_path, motif_length, num_cores):
    """Combine all *_motifs.tsv files into a motifs x libraries DataFrame."""
    input_dir = Path(directory_path)
    file_paths = sorted(input_dir.glob("*_motifs.tsv"))

    if not file_paths:
        raise FileNotFoundError(
            f"No *_motifs.tsv files found in {input_dir}"
        )

    motifs = all_motifs(motif_length)
    motifs_set = set(motifs)
    libraries = [sample_name(path) for path in file_paths]

    motif_df = pd.DataFrame(0, index=motifs, columns=libraries, dtype="int64")
    motif_df.index.name = "motif"

    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        future_to_path = {
            executor.submit(read_motif_counts, path, motifs_set): path
            for path in file_paths
        }

        for future in tqdm(
            as_completed(future_to_path),
            total=len(future_to_path),
            desc="Processing motif count files",
        ):
            path = future_to_path[future]
            library = sample_name(path)
            counts = future.result()

            for motif, count in counts.items():
                motif_df.at[motif, library] = count

    return motif_df


def main():
    if len(sys.argv) != 5:
        print(
            "Usage: python motif_frequency_counter_efficient.py "
            "<directory_of_sample_files> <motif_length> <num_cores> "
            "<output_directory>",
            file=sys.stderr,
        )
        sys.exit(1)

    directory_path = sys.argv[1]
    motif_length = int(sys.argv[2])
    num_cores = int(sys.argv[3])
    output_directory = Path(sys.argv[4])

    output_directory.mkdir(parents=True, exist_ok=True)

    start_time = time.time()
    motif_counts_df = combine_motif_counts(
        directory_path,
        motif_length,
        num_cores,
    )

    output_counts_file = output_directory / "motif_counts.tsv"
    output_relfreq_file = output_directory / "motifs_rel_freq_wide.tsv"

    # Counts: rows = motifs, columns = libraries.
    motif_counts_df.to_csv(output_counts_file, sep="\t")

    # Relative frequencies: rows = libraries, columns = motifs.
    # OTHER is excluded, so each non-empty library is normalized over the
    # canonical A/C/G/T motifs only.
    column_sums = motif_counts_df.sum(axis=0).replace(0, pd.NA)
    relative_frequency_df = motif_counts_df.divide(column_sums, axis=1).T
    relative_frequency_df.index.name = "library"
    relative_frequency_df = relative_frequency_df.sort_index()
    relative_frequency_df.to_csv(output_relfreq_file, sep="\t")

    elapsed = time.time() - start_time
    print(f"Processed in {elapsed:.2f} seconds.")
    print(f"Motif counts saved to {output_counts_file}")
    print(f"Relative frequencies saved to {output_relfreq_file}")


if __name__ == "__main__":
    main()
