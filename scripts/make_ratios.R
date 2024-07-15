#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
# For single cutpoint pathway, input_file looks like "frag_counts_by_len_class.tsv"
# For multiple cutpoint pathway, input_file looks like "frag_counts_{cutpoint}.tsv"
frags_tsv <- args[1] # fragment counts file with all short/long fragment lengths for all libraries
ratios_tsv <- args[2] # output file
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

# Load aggregate frag tsv
cat("Loading frag file...\n")
frags <- read_tsv(frags_tsv)

# From per-position, per library short and long fragment counts, zero-centered fragment ratio
# See https://github.com/cancer-genomics/reproduce_lucas_wflow/blob/master/analysis/fig2a.Rmd

cat("Calculating short to long fragment ratios...\n")
ratios <-
  frags %>%
  mutate_at(vars(start, end, count), as.numeric) %>%
  # Put lib-bin short and long values on same row in order to make per-row ratios
  pivot_wider(names_from = len_class, values_from = count) %>%
  mutate(fract = short / long) %>%
  select(library, chr, start, end, fract) %>%
  # Zero center by library
  group_by(library) %>%
  mutate(ratio.centered = scale(fract, center = TRUE, scale = FALSE)[, 1])

cat("Writing ratios file...\n")
write_tsv(ratios, file = ratios_tsv)

cat("Script completed successfully")