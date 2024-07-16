#!/usr/bin/env Rscript

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
# args[1] is a space-separated list of fragment count file names. 
# counts_files is a vector of file names.
counts_files <- unlist(strsplit(args[1], " ")) 
frag_counts <- args[2] #output
num_threads <- as.numeric(args[3]) #parameter
log_file <- args[4]

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

# Function to read and process total files
process_total_files <- function(files) {
  cat("Processing total files: ", length(files), " files...\n")
  dt_list <- lapply(files, function(file) {
    dt <- fread(file, header = FALSE, col.names = c("chr", "start", "end", "gc", "count"),
                colClasses = c("character", "integer", "integer", "numeric", "integer"))
    dt[, library := sub("_count_total.tsv", "", basename(file))]
    return(dt)
  })
  dt_total <- rbindlist(dt_list)
  setorder(dt_total, library, chr, start, end)
  dt_total <- dt_total[, .(library, chr, start, end, gc, count)]
  fwrite(dt_total, frag_counts, sep = "\t", col.names = TRUE)
}

# Process total files
process_total_files(counts_files)

cat("Script completed successfully\n")