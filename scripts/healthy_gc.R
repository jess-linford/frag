#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
healthy_libs_str <- args[1] # Space-separated list of healthy library gc_distro file names
healthy_med_file <- args[2] # output
log_file <- args[3]

library(tidyverse)

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

cat("Creating list of healthy libraries...\n")
healthy_libs_distros <- unlist(strsplit(healthy_libs_str, " "))

read_in_gc <- function(gc_csv){
  read.csv(gc_csv, header = TRUE)
}

cat("Reading in gc distributions of healthy libraries...\n")
healthy_list <- lapply(healthy_libs_distros, read_in_gc)

# Bind
healthy_all <- do.call(rbind, healthy_list)

# Calculate library sizes
cat("Calculating target library size...\n")
lib_sizes <- healthy_all %>%
  group_by(library_id) %>%
  summarise(lib_size = sum(frag_count))

# Calculate target_lib_size = median(library_size) for healthy distributions
target_lib_size = median(lib_sizes$lib_size)

# Create target distribution:
# Target frag count in each stratum = median(frag fraction) * median(lib size)
# Not sure if this is needlessly complicated - better to just take median(frag count)?
# Results are similar, but not identical
cat("Creating target median gc distribution...\n")
healthy_med <-
  healthy_all %>%
  group_by(gc_strata) %>%
  summarise(med_frag_fract = median(frag_fract)) %>%
  mutate(target_frag_count = round(med_frag_fract * target_lib_size))

# Save
cat("Writing output file...\n")
saveRDS(healthy_med, file = healthy_med_file)

cat("Script completed successfully")