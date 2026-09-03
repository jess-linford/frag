#!/usr/bin/env bash

set -euo pipefail

input=$1
threads=$2
output=$3

# Filter to reads that are:
#  - Paired with both read and mate mapped
#  - Primary alignments only
#  - Non-duplicate
#  - MAPQ >= 20
#
# DO NOT require the "proper pair" flag, because this may
# remove valid cfDNA fragments.

samtools view \
    -@ "$threads" \
    -b \
    -F 3340 \
    -q 20 \
    -o "$output" \
    "$input"

samtools index -@ "$threads" "$output"