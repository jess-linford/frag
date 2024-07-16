#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
# args[1] and args[2] are space-separated lists of short and long fragment count file names.
# short_input_files and long_input_files are vectors of file names.
short_input_files <- unlist(strsplit(args[1], " "))
long_input_files <- unlist(strsplit(args[2], " "))
# args[3] is a space-separated list of frag_counts_{cutpoint}.tsv output file names.
# output_files is a vector of file names.
output_files <- unlist(strsplit(args[3], " "))
num_threads <- as.integer(args[4])
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

# Function to read, process, and sort counts from files
process_counts <- function(files, len_class) {
  counts_list <- lapply(files, function(file) {
    dt <- fread(file, header = FALSE, col.names = c("chr", "start", "end", "gc", "count"),
                colClasses = c("character", "integer", "integer", "numeric", "integer"))
    dt[, library := sub(paste0("_count_", len_class, "_\\d+\\.tsv$"), "", basename(file))]
    dt[, len_class := len_class]
    return(dt)
  })
  combined_counts <- rbindlist(counts_list)
  return(combined_counts)
}

# Loop over output files to combine counts
for (i in seq_along(output_files)) {
  cutpoint <- sub("frag_counts_(\\d+)\\.tsv", "\\1", basename(output_files[i]))
  
  short_files <- grep(paste0("_count_short_", cutpoint, "\\.tsv$"), short_input_files, value = TRUE)
  long_files <- grep(paste0("_count_long_", cutpoint, "\\.tsv$"), long_input_files, value = TRUE)
  
  cat(sprintf("Processing counts for cutpoint %s...\n", cutpoint))
  short_counts <- process_counts(short_files, "short")
  long_counts <- process_counts(long_files, "long")
  
  combined_counts <- rbind(short_counts, long_counts)
  
  # Sort rows and select columns before writing to file
  setorder(combined_counts, library, len_class, chr, start, end)
  combined_counts <- combined_counts[, .(library, len_class, chr, start, end, gc, count)]
  
  cat(sprintf("Writing combined counts to %s...\n", output_files[i]))
  fwrite(combined_counts, file = output_files[i], sep = "\t", col.names = TRUE)
}

cat("Script completed successfully\n")