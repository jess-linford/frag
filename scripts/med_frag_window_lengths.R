#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
frag_bed_file <- args[1] #input
keep_bed_file <- args[2] #input
frag_length_low <- as.numeric(args[3]) #parameter
frag_length_high <- as.numeric(args[4]) #parameter
num_threads <- as.integer(args[5]) #parameter
med_window_lengths <- args[6] #output
log_file <- args[7]

library(data.table)
library(GenomicRanges)
setDTthreads(num_threads)

# Open a sink connection to capture all output
log_conn <- file(log_file, open = "wt")
sink(log_conn)
sink(log_conn, type = "message")

# Ensure sinks are closed upon exit
on.exit({
  sink()
  sink(type = "message")
  close(log_conn)
})

# Read and filter the fragment bed file
cat("Loading input fragment bed file...\n")
frag_bed <- fread(frag_bed_file,
                  col.names = c("chr", "start", "end", "gc_raw", "len"),
                  colClasses = c("character", "integer", "integer", "numeric", "integer"))

cat("Filtering to short autosomal reads...\n")
frag_bed <- frag_bed[len >= frag_length_low & len <= frag_length_high & chr %in% paste0("chr", seq(1:22))]

setorder(frag_bed, chr, start, end)

frag_bed <- frag_bed[, .(chr, start, end, len)]

# Read the keep bed file
cat("Reading 5Mb keep bed file...\n")
keep_bed <- fread(keep_bed_file, col.names = c("chr", "start", "end", "gc"))

# Create GenomicRanges objects for overlap calculation
cat("Creating GenomicRanges objects...\n")
input_gr <- GRanges(seqnames = frag_bed$chr,
                    ranges = IRanges(start = frag_bed$start, end = frag_bed$end))
keep_bed_gr <- GRanges(seqnames = keep_bed$chr,
                       ranges = IRanges(start = keep_bed$start, end = keep_bed$end))

cat("Finding overlaps...\n")
hits <- findOverlaps(keep_bed_gr, input_gr)

frag_bed <- frag_bed[subjectHits(hits), window_id := queryHits(hits)]

# Calculate median lengths
cat("Calculating median lengths...\n")
median_lengths <- frag_bed[, .(med_len = median(len)), by = window_id]

# Merge with keep_bed to get the final result
cat("Creating median window lengths table...\n")
keep_bed$window_id <- 1:nrow(keep_bed)
median_lengths_dt <- merge(keep_bed, median_lengths, by = "window_id", all.x = TRUE)

cat("Writing median window lengths file...\n")
fwrite(median_lengths_dt[, .(chr, start, end, med_len)], file = med_window_lengths, sep = "\t", col.names = FALSE)

cat("Script completed successfully\n")