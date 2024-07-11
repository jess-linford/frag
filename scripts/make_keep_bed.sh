# Assign input arguments to variables
gc5mb=$1
blklist=$2
keep=$3

# Filter gc5mb file for desired chromosomes, intersect with blklist, and keep rows with gc >= 0.3
awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/' "$gc5mb" | \
bedtools intersect -a - -b "$blklist" -v -wa | \
grep -v _ | \
awk '{ if ($4 >= 0.3 && $4 != "nan") print $0 }' | \
bedtools sort -i - > "$keep"
