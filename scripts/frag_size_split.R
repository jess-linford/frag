#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
frag_weights_bed <- args[1] #input
cutpoint <- as.numeric(args[2]) #parameter
short_bed_file <- args[3] #output
long_bed_file <- args[4] #output
num_threads <- args[5] #parameter
log_file <- args[6]

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

# Load the fragment weights bed file
cat("Loading fragment weights bed file...\n")
frag <- fread(frag_weights_bed,
              col.names = c("chr", "start", "end", "len", "gc_strata", "weight"),
              colClasses = c("character", "integer", "integer", "integer", "numeric", "numeric"))

# Separate short and long fragments
cat("Creating short and long fragment datasets...\n")
frag_short <- frag[len <= cutpoint]
frag_long <- frag[len > cutpoint]

# Write short and long fragment bedfiles
cat("Writing short and long fragment bedfiles...\n")
fwrite(frag_short, file = short_bed_file, sep = "\t", col.names = FALSE)
fwrite(frag_long, file = long_bed_file, sep = "\t", col.names = FALSE)

cat("Script completed successfully\n")