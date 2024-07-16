#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
healthy_med_gc_distro <- args[1] #input
frag_bed_file <- args[2] #input
num_threads <- as.integer(args[3]) #parameter
frag_length_low <- as.numeric(args[4]) #parameter
cutpoints <- as.numeric(unlist(strsplit(args[5], ","))) #parameter
frag_length_high <- as.numeric(args[6]) #parameter
frag_weights_bed <- args[7] #output
log_file <- args[8]

cat("Parsed cutpoints: ", paste(cutpoints, collapse = ", "), "\n")

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

# Load the healthy median gc distribution
cat("Loading healthy median gc distribution...\n")
healthy_med <- readRDS(healthy_med_gc_distro)

# Load the input fragment bed file
cat("Loading input fragment bed file...\n")
frag_bed <- fread(frag_bed_file,
                  col.names = c("chr", "start", "end", "gc_raw", "len"),
                  colClasses = c("character", "integer", "integer", "numeric", "integer"))

# Calculate GC strata and filter to reads on autosomes with lengths in given range
cat("Calculating GC strata and filtering to short autosomal reads...\n")
frag_bed[, gc_strata := round(gc_raw, 2)]
frag_bed <- frag_bed[len >= frag_length_low & len <= frag_length_high & chr %in% paste0("chr", seq(1:22))]

# Join fragment bed file to healthy median gc distribution
cat("Merging bed file to healthy median gc distribution...\n")
frag <- merge(frag_bed, healthy_med, by = "gc_strata", all.x = TRUE)

# Calculate weight = target_frag_count / num_frags_in_stratum
cat("Calculating fragment weights...\n")
frag[, n := .N, by = gc_strata]
frag[, weight := target_frag_count / n]

# Sort the sampled data.table by the start and end columns
setorder(frag, chr, start, end)

# Select the required columns
frag <- frag[, .(chr, start, end, len, gc_strata, weight)]

# Write fragment weights bedfile
cat("Writing fragment weights bedfile...\n")
fwrite(frag, file = frag_weights_bed, sep = "\t", col.names = FALSE)

base_name <- sub("_frag.bed$", "", basename(frag_bed_file))

# Loop through each cutpoint and create short and long fragment datasets
for (cutpoint in cutpoints) {
  cat(sprintf("Processing cutpoint: %d\n", cutpoint))
  short_frag_bed <- paste0(dirname(frag_bed_file), "/", base_name, "_short_weights_", cutpoint, ".bed")
  long_frag_bed <- paste0(dirname(frag_bed_file), "/", base_name, "_long_weights_", cutpoint, ".bed")
  
  cat(sprintf("Creating short and long fragment datasets for cutpoint %d...\n", cutpoint))
  frag_short <- frag[len <= cutpoint]
  frag_long <- frag[len > cutpoint]

  cat(sprintf("Short fragment dataset size: %d\n", nrow(frag_short)))
  cat(sprintf("Long fragment dataset size: %d\n", nrow(frag_long)))
  
  cat("Writing file:", short_frag_bed, "...\n")
  fwrite(frag_short, file = short_frag_bed, sep = "\t", col.names = FALSE)
  cat("Writing file:", long_frag_bed, "...\n")
  fwrite(frag_long, file = long_frag_bed, sep = "\t", col.names = FALSE)
}

cat("Script completed successfully\n")