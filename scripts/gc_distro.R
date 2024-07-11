#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
bed_file <- args[1] #input
distro_file <- args[2] #output
frag_length_low <- as.integer(args[3])
frag_length_high <- as.integer(args[4])
num_threads <- as.integer(args[5])
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

cat("Arguments:\n")
cat("bed_file:", bed_file, "\n")
cat("distro_file:", distro_file, "\n")
cat("frag_length_low:", frag_length_low, "\n")
cat("frag_length_high:", frag_length_high, "\n")
cat("num_threads:", num_threads, "\n")
cat("log_file:", log_file, "\n")

cat("Reading input BED file:", bed_file, "\n")
bed <- fread(bed_file,
             col.names = c("chr", "start", "end", "gc_raw", "len"),
             colClasses = c("character", "integer", "integer", "numeric", "integer"))

cat("BED file read successfully. Number of rows: ", nrow(bed), "\n")

# Generate distribution csv
cat("Generating distribution dataframe...\n")
# Create rounded gc column
bed[, gc_strata := round(gc_raw, 2)]
# Filter the bed file to reads on autosomes with lengths in given range
bed <- bed[len >= frag_length_low & len <= frag_length_high & chr %in% paste0("chr", seq(1:22))]
cat("Number of rows after filtering:", nrow(bed), "\n")

cat("Creating distribution dataframe...\n")
# Count number of rows for each gc_strata value and store result as 'frag_count'
distro <- bed[, .(frag_count = .N), by = gc_strata]
# Calculate proportion of fragments in each gc_strata value
distro[, frag_fract := frag_count / sum(frag_count)]
# Add library id
distro[, library_id := gsub("_frag\\.bed$", "", basename(bed_file))]

# Sort the data by library_id and gc_strata
setorder(distro, library_id, gc_strata)
cat("Number of rows in distribution dataframe:", nrow(distro), "\n")

cat("Writing distribution file...\n")
# Select and write the required columns
fwrite(distro[, .(library_id, gc_strata, frag_count, frag_fract)], file = distro_file, row.names = FALSE)

cat("Script completed successfully\n")
