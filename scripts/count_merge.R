#!/usr/bin/env Rscript

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
# args[1] is a space-separated list of fragment count file names. 
# counts_files is a vector of file names.
counts_files <- unlist(strsplit(args[1], " ")) 
frag_counts <- args[2] #output
frag_counts_by_len_class <- args[3] #output
num_threads <- as.numeric(args[4]) #parameter
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

# Function to read and process short and long files
process_len_class_files <- function(files, len_class) {
  cat("Processing ", len_class, " files: ", length(files), " files...\n", sep = "")
  dt_list <- lapply(files, function(file) {
    dt <- fread(file, header = FALSE, col.names = c("chr", "start", "end", "gc", "count"),
                colClasses = c("character", "integer", "integer", "numeric", "integer"))
    dt[, library := sub(paste0("_count_", len_class, ".tsv"), "", basename(file))]
    dt[, len_class := len_class]
    return(dt)
  })
  return(dt_list)
}

# Separate the input files based on their patterns
total_files <- grep("total", counts_files, value = TRUE)
short_files <- grep("short", counts_files, value = TRUE)
long_files <- grep("long", counts_files, value = TRUE)

# Process total files
process_total_files(total_files)

# Process short and long files
dt_short <- process_len_class_files(short_files, "short")
dt_long <- process_len_class_files(long_files, "long")

# Concatenate and order the len_class files
dt_len_class <- rbindlist(c(dt_short, dt_long))
setorder(dt_len_class, library, len_class, chr, start, end)
dt_len_class <- dt_len_class[, .(library, len_class, chr, start, end, gc, count)]
fwrite(dt_len_class, frag_counts_by_len_class, sep = "\t", col.names = TRUE)

cat("Script completed successfully\n")