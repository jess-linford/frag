#!/usr/bin/env bash

# Snakemake variables
outdir="$1"
input_bam="$2"
output_metrics="$3"
output_histogram="$4"
max_width="$5"

# Create the output directory
mkdir -p "$outdir"

# Run the picard CollectInsertSizeMetrics command
picard CollectInsertSizeMetrics \
-I "$input_bam" \
-O "${output_metrics}.tmp" \
-H "$output_histogram" \
-W "$max_width"

# Remove the header rows and keep only the table with sizes and counts
awk 'NR > 11 {print}' "${output_metrics}.tmp" > "$output_metrics"
rm "${output_metrics}.tmp"