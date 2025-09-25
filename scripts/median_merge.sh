#!/usr/bin/env bash
set -euo pipefail

medians_dir="${1}"
out_long="${2}"
out_wide="${3}"

# Remove old outputs if present
[ -f "$out_long" ] && rm "$out_long"
[ -f "$out_wide" ] && rm "$out_wide"

# Build the long table
for file in "${medians_dir}"/*_med_frag_window_lengths.tsv; do
    [ -e "$file" ] || { echo "No files matched in ${medians_dir}"; exit 1; }
    library=$(basename "$file" "_med_frag_window_lengths.tsv")
    awk -v lib="$library" 'BEGIN {OFS="\t"} {print lib, $0}' "$file" >> "$out_long"
done

# Add a header to the long table
sed -i '1 i\library\tchr\tstart\tend\tmedian' "$out_long"

# Transform to wide using data.table
Rscript - <<RS
suppressMessages(library(data.table))
dt <- fread("$out_long")
wide <- dcast(dt, library ~ chr + start + end, value.var = "median")
setorder(wide, library)
fwrite(wide, file = "$out_wide", sep = "\t", quote = FALSE, na = "NA")
RS

echo "Wrote long table to: $out_long"
echo "Wrote wide table to: $out_wide"
