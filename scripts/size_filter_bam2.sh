#!/usr/bin/env bash

set -euo pipefail

input=$1
threads=$2
min_length=$3
max_length=$4
output=$5

echo "========================================"
echo "Fragment length BAM filtering"
echo "========================================"
echo "Date: $(date)"
echo "Hostname: $(hostname)"
echo
echo "Input BAM:     $input"
echo "Output BAM:    $output"
echo "Threads:       $threads"
echo "Min length:    $min_length"
echo "Max length:    $max_length"
echo

echo "Samtools:"
which samtools
samtools --version
echo

echo "Checking input BAM..."
ls -lh "$input"
samtools quickcheck -v "$input"

echo
echo "Creating output directory..."
mkdir -p "$(dirname "$output")"

echo
echo "Running fragment-length filter..."

samtools view \
    -@ "$threads" \
    -f 2 \
    -F 3328 \
    -e "((tlen >= $min_length && tlen <= $max_length) || (tlen <= -$min_length && tlen >= -$max_length))" \
    -O BAM \
    -o "$output" \
    "$input"

echo
echo "Filtering complete."

echo "Checking output BAM..."
ls -lh "$output"
samtools quickcheck -v "$output"

echo
echo "Indexing output BAM..."

samtools index \
    -@ "$threads" \
    "$output"

echo
echo "Checking index..."
ls -lh "${output}.bai"

echo
echo "Number of output alignment records:"
samtools view -c "$output"

echo
echo "Fragment-length range in output:"
samtools view "$output" |
    awk '$9 > 0 {
        if (min == "" || $9 < min) min=$9
        if (max == "" || $9 > max) max=$9
        n++
    }
    END {
        print "Positive-TLEN fragments:", n+0
        print "Minimum TLEN:", min
        print "Maximum TLEN:", max
    }'

echo
echo "Completed successfully."