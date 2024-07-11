#!/usr/bin/env bash

# Snakemake variables
input_bam="$1"
output="$2"
max_width="$3"

# Run the picard CollectInsertSizeMetrics command
picard CollectInsertSizeMetrics \
-I "$input_bam" \
-O "${output}.tmp" \
-W "$max_width"

# Remove the header rows and keep only the table with sizes and counts
awk 'NR > 11 {print}' "${output}.tmp" > "$output"
rm "${output}.tmp"