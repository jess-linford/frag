#!/usr/bin/env bash

medians_dir="${1}"
out_tsv="${2}"

# Remove the existing aggregate file if present
if [ -f "$out_tsv" ]; then 
    rm "$out_tsv"
fi

# Make aggregate file
for file in "${medians_dir}"/*_med_frag_window_lengths.tsv; do
    library=$(basename "$file" "_med_frag_window_lengths.tsv")
    awk -v lib="$library" 'BEGIN {OFS="\t"} {print lib, $0}' "$file" >> "$out_tsv"
done

# Add a header
sed -i '1 i\library\tchr\tstart\tend\tmedian' "$out_tsv"