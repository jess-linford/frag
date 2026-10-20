#!/usr/bin/env bash

# Extract strand-aware 5' fragment-end motifs from a paired-end BAM.
#
# For each selected alignment:
#   + strand: take N reference bases starting at the alignment start.
#   - strand: take N reference bases ending at the alignment end.
#
# bedtools getfasta -s respects the BED strand, so sequence from a minus-strand
# interval is reverse-complemented. This reports every motif in the biological
# 5' -> 3' orientation of the aligned read rather than reference orientation.
#
# Usage:
#   extract_5p_motifs.sh <bam> <fasta> <out_counts.tsv> <motif_length> [threads] [region]
#
# Output:
#   motif<TAB>count
#
# An OTHER row contains extracted sequences that are not exactly motif_length
# bases of A/C/G/T (for example, sequences containing N).
#
# Read filtering:
#   MAPQ >= 30
#   paired
#   mapped
#   primary alignment
#   non-duplicate
#
# Both R1 and R2 are included. Alignment strand, not mate identity, determines
# which genomic edge is treated as the biological 5' end.
#
# samtools and bedtools are expected to be available on PATH (for example from
# an activated Conda environment).

set -euo pipefail

bam="$1"
fasta="$2"
out="$3"
motif_length="$4"
threads="${5:-3}"
region="${6:-}"

# motif_length must be a positive integer because it is used to construct BED
# coordinates and to validate the extracted sequence length.
if ! [[ "$motif_length" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: motif_length must be a positive integer; got '$motif_length'." >&2
    exit 1
fi

# Build the optional samtools region argument as an array so that an empty
# region does not create a spurious argument.
region_args=()
if [[ -n "$region" ]]; then
    region_args=("$region")
fi

samtools view \
    -b \
    -q 30 \
    -f 1 \
    -F 3340 \
    --threads "$threads" \
    "$bam" \
    "${region_args[@]}" |
bedtools bamtobed -i stdin |
awk -v n="$motif_length" '
    BEGIN {
        OFS = "\t"
    }

    # BED columns from bedtools bamtobed:
    #   $1 = chromosome
    #   $2 = start (0-based)
    #   $3 = end (exclusive)
    #   $6 = alignment strand

    # Plus strand: the 5-prime end is the lower genomic coordinate.
    # Keep [start, start + n).
    $6 == "+" {
        print $1, $2, $2 + n, ".", ".", $6
    }

    # Minus strand: the 5-prime end is the higher genomic coordinate.
    # Keep [end - n, end). The coordinate check prevents a negative start.
    $6 == "-" {
        if ($3 >= n) {
            print $1, $3 - n, $3, ".", ".", $6
        }
    }
' |
bedtools getfasta \
    -fi "$fasta" \
    -bed stdin \
    -s \
    -tab |
awk -F '\t' -v n="$motif_length" '
    {
        motif = toupper($2)

        # Count only motifs with exactly n canonical DNA bases. Avoid an
        # interval-regex such as {n} here for compatibility across awk variants.
        if (length(motif) == n && motif !~ /[^ACGT]/) {
            count[motif]++
        } else {
            other++
        }
    }

    END {
        for (motif in count) {
            print motif "\t" count[motif]
        }
        print "OTHER\t" other + 0
    }
' |
sort > "$out"
