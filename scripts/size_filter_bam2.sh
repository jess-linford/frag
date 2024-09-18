#!/usr/bin/env bash

input=$1
threads=$2
min_length=$3
max_length=$4
output=$5

# Filter BAM file for reads between min_length and max_length and write to output file
samtools view -@ $threads -e "length(seq) >= $min_length && length(seq) <= $max_length" -O BAM -o $output $input

# Index the final output
samtools index $output