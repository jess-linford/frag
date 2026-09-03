#!/usr/bin/env bash

set -euo pipefail

input=$1
threads=$2
max_length=$3
output=$4

# Filter to reads that are:
# - Paired, with both read and mate mapped
# - Primary alignments only
# - Non-duplicate
# - MAPQ >= 20
# - Fragment/template length <= max_length
# - On autosomes chr1-chr22
#
# DO NOT require the "proper pair" flag.

regions=$(seq -f "chr%.0f" 1 22)

samtools view \
    -@ "$threads" \
    -b \
    -f 1 \
    -F 3340 \
    -q 20 \
    -e "((tlen > 0 && tlen <= $max_length) || (tlen < 0 && tlen >= -$max_length))" \
    -o "$output" \
    "$input" \
    $regions

samtools index -@ "$threads" "$output"