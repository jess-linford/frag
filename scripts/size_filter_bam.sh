#!/usr/bin/env bash

input=$1
threads=$2
max_length=$3
output=$4

# Filter BAM file for reads <= 500bp and write to output file
samtools view -@ $threads -e "length(seq) <= $max_length" -O BAM -o $output $input

# Index the final output
samtools index $output