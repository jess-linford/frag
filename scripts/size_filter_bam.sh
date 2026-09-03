#!/usr/bin/env bash

set -euo pipefail

input=$1
threads=$2
max_length=$3
output=$4

# Filter BAM to fragments with inferred template length <= max_length.
# TLEN has opposite signs for the two mates, so filter on |TLEN|.
# Exclude TLEN = 0 because fragment length is undefined/uninformative there.

samtools view \
    -@ "$threads" \
    -e "((tlen > 0 && tlen <= $max_length) || (tlen < 0 && tlen >= -$max_length))" \
    -O BAM \
    -o "$output" \
    "$input"

samtools index -@ "$threads" "$output"