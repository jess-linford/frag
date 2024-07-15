#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
# For single cutpoint pathway, input_file looks like {library}_{short/long}_weights.bed
# or {library}_weights.bed
# For multiple cutpoints pathway, input_file looks like 
# {library}_{short/long}_weights_{cutpoint}.bed
input_file <- args[1]
keep_bed_file <- args[2] # keep_5mb.bed
# For single cutpoint pathway, output_file looks like {library}_count_{short/long/total}.tmp
# For multiple cutpoints pathway, output_file looks like {library}_count_{short/long}_{cutpoint}.tmp
output_file <- args[3]
num_threads <- as.integer(args[4])
log_file <- args[5]

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

library(data.table)
library(GenomicRanges)
setDTthreads(num_threads)

# Read input files with specified column names
cat("Reading fragment weight bed file...\n")
input <- fread(input_file, col.names = c("chr", "start", "end", "len", "gc_strata", "weight"))
cat("Reading 5mb keep bed file...\n")
keep_bed <- fread(keep_bed_file, col.names = c("chr", "start", "end", "gc"))

# Sort the input file
setorder(input, chr, start, end)

# Create GenomicRanges objects
cat("Creating GenomicRanges objects...\n")
input_gr <- GRanges(seqnames = input$chr,
                    ranges = IRanges(start = input$start, end = input$end),
                    weight = input$weight)
keep_bed_gr <- GRanges(seqnames = keep_bed$chr,
                       ranges = IRanges(start = keep_bed$start, end = keep_bed$end))

# Use GenomicRanges findOverlaps and sum weights
cat("Finding overlaps...\n")
hits <- findOverlaps(keep_bed_gr, input_gr)
cat("Counting fragments by summing weights...\n")
keep_bed$count <- tapply(input$weight[subjectHits(hits)], queryHits(hits), sum, default = 0)

# Create the output data.table
cat("Creating output...\n")
output_dt <- data.table(chr = as.character(seqnames(keep_bed_gr)),
                        start = start(keep_bed_gr),
                        end = end(keep_bed_gr),
                        gc = keep_bed$gc,
                        count = as.integer(round(keep_bed$count)))

# Write the output file
cat("Writing output file...\n")
fwrite(output_dt, file = output_file, sep = "\t", col.names = FALSE)

cat("Script completed successfully")
