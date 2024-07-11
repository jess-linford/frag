#!/usr/bin/env bash

input=$1
threads=$2
max_length=$3
output=$4

# Filter to reads that are
#  - Excluding any unmapped, not primary alignment, or duplicates
#  - Only MAPQ > 20
#  - <= max_length
#  - On autosomes
# DO NOT restrict to "proper pairs"- this clips long cfDNA fragments!

# Define the regions for autosomes (chromosomes 1-22)
regions=$(seq -f "chr%.0f" 1 22)

# Filter the bam file
samtools view -@ $threads -b -F 1284 -q 20 -e "length(seq) <= $max_length" -o $output $input $regions

# Index the final output
samtools index $output