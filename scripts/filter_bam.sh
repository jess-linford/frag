#!/usr/bin/env bash

set -euo pipefail

input=$1
threads=$2
max_length=$3
output=$4

# Retain paired, primary, non-duplicate alignments where:
# - both reads are mapped
# - MAPQ >= 20
# - inferred fragment/template length <= max_length
#
# Do NOT require the SAM proper-pair flag, because that may
# remove biologically valid longer cfDNA fragments.

samtools view \
    -@ "$threads" \
    -b \
    -f 1 \
    -F 3340 \
    -q 20 \
    -e "((tlen > 0 && tlen <= $max_length) || (tlen < 0 && tlen >= -$max_length))" \
    -o "$output" \
    "$input"

samtools index -@ "$threads" "$output"