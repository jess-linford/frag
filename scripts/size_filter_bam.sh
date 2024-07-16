#!/usr/bin/env bash

input=$1
threads=$2
output=$3

# Filter BAM file for reads <= 500bp and write to output file
samtools view -@ $threads -e 'length(seq) <= 500' -O BAM -o $output $input

# Index the final output
samtools index $output