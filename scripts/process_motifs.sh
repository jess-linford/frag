#!/bin/bash

in_bam="${1}"
in_fasta="${2}"
n_motif="${3}"
threads="${4}"
out_merged="${5}"

# Function to create index if not present
check_and_index_bam() {
    local bam_file="${1}"
    if [ ! -f "${bam_file}.bai" ]; then
        echo "Index file for ${bam_file} not found. Creating index..."
        samtools index "${bam_file}"
        if [ $? -ne 0 ]; then
            echo "Failed to create index for ${bam_file}. Exiting."
            exit 1
        fi
    fi
}

main(){
    # Debugging: Print all inputs
    echo "Input BAM: $in_bam"
    echo "FASTA: $in_fasta"
    echo "Motif length: $n_motif"
    echo "Threads: $threads"
    echo "Output file: $out_merged"

    # Check for missing inputs
    if [ -z "$in_bam" ] || [ -z "$in_fasta" ] || [ -z "$n_motif" ] || [ -z "$threads" ] || [ -z "$out_merged" ]; then
        echo "Error: One or more input parameters are missing."
        exit 1
    fi

    # Check and create BAM index if needed
    check_and_index_bam "$in_bam"

    # Extract motifs
    forward_motif \
        "$in_bam" \
        "$threads" \
        "$in_fasta" \
        "$n_motif" > "$out_merged"

    reverse_motif \
        "$in_bam" \
        "$threads" \
        "$in_fasta" \
        "$n_motif" >> "$out_merged"
}

forward_motif(){
    #
    local in_bam="${1}"
    local threads="${2}"
    local in_fasta="${3}"
    local n_motif="${4}"
    #
    # Take first read in mapped, paired, with normal FS orientation.
    # View perfect matching reads (for BWA), first in pair.
    samtools view \
             -h \
             -q 30 \
             -f 65 \
             --threads "$threads" "$in_bam" |
        # Fetch reference
        bedtools bamtobed -i stdin |
        bedtools getfasta -bed stdin -fi "$in_fasta" |
        # Sed magic to extract motifs from fasta
        sed "1d; n; d" | sed -E "s/(.{$n_motif}).*/\1/"
}

reverse_motif(){
    #
    local in_bam="${1}"
    local threads="${2}"
    local in_fasta="${3}"
    local n_motif="${4}"
    #
    # Take SECOND read in mapped, paired, with normal FS orientation.
    # View perfect matching reads (for BWA).
    samtools view \
             -h \
             -q 30 \
             -f 129 \
             --threads "$threads" "$in_bam" |
        # Fetch reference
        bedtools bamtobed -i stdin |
        bedtools getfasta -bed stdin -fi "$in_fasta" |
        # Sed magic to extract motifs from fasta
        sed "1d; n; d" | sed -E "s/.*(.{$n_motif})/\1/" |
        # Generate reverse complement
        tr ACGT TGCA | rev
}

main "$@"
