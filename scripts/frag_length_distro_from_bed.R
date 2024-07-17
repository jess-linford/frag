#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
bed_file <- args[1] #input bed file
length_distro_file <- args[2] #output file
max_frag_length <- as.integer(args[3]) #parameter
num_threads <- as.integer(args[4]) #parameter
log_file <- args[5]

library(data.table)
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

# Load the input fragment bed file
cat("Reading fragment bed file...\n")
bed <- fread(bed_file, col.names = c("chr", "start", "end", "gc_raw", "len"),
             colClasses = c("character", "integer", "integer", "numeric", "integer"))

bed <- bed[len <= max_frag_length]

# Create a table with fragment length and count of reads
cat("Creating fragment length distribution table...\n")
length_table <- bed[, .(count = .N), by = .(len)]
setorder(length_table, len, count)

# Select and write the required columns
cat("Writing fragment length distribution table to file...\n")
fwrite(length_table[, .(len, count)], file = length_distro_file, row.names = FALSE)

cat("Script completed successfully\n")