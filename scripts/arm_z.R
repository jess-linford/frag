#########1#########2#########3#########4#########5#########6#########7#########8
# Mathios, et al. 2021. https://doi.org/10.1038/s41467-021-24994-w
# "In addition to fragmentation profiles, we also computed z-scores for
#  chromosomal arms and a genome-wide summary of the overall cfDNA fragmentation
#  as previously described23,51 with the modifications indicated below. We only
#  analyzed the 39 chromosomal arms that were not acrocentric. Z-scores for each
#  of the 39 autosomal arms were obtained by centering and scaling the total
#  GC-adjusted fragment count for each arm by the mean and standard deviation of
#  the corresponding arm-specific counts in the 54 non-cancer samples used as a
#  reference set."

args <- commandArgs(trailingOnly = TRUE)
cytobands_tsv <- args[1] #input
libraries_file <- args[2] #input - must have columns "library" and "cohort"
frag_counts_tsv <- args[3] #input
armz_tsv <- args[4] #output
log_file <- args[5]

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

# Read input files
cat("Reading in cytobands file...\n")
cytobands <- read_tsv(cytobands_tsv, col_names = c("chr", "start", "end", "band", "stain"))
cat("Reading in frag_counts file...\n")
frag_counts <- read_tsv(frag_counts_tsv)
cat("Reading in libraries file...\n")
libraries <- read_tsv(libraries_file)

# If you don't care about duplicate healthy libraries, comment out line 44
cat("Getting list of healthy cohort libraries...\n")
healthy_libs <- libraries %>% filter(cohort == "healthy") %>% 
  # filter(duplicate == 0 | (duplicate == 1 & batch != 1)) %>%
  pull(library)

# Define non-acrocentric arms and their genomic coordinates
cat("Creating dataframe of non-acrocentric arms with coordinates...\n")
autosomes <- c(paste0("chr", seq(1:22)))
noacro_chroms <- c(paste0("chr", seq(1:12)), paste0("chr", seq(16, 20, by = 1)))
acro_chroms <- autosomes[!autosomes %in% noacro_chroms]
noacro <- cytobands %>%
  mutate(arm = substr(band, 1, 1)) %>%
  group_by(chr, arm) %>%
  summarise(armstart = min(start),
            armend = max(end)) %>%
  filter(chr %in% noacro_chroms | (chr %in% acro_chroms & arm == "q")) %>%
  pivot_wider(names_from = arm, values_from = c(armstart, armend))

# Find total number of counts by library
counts_by_lib <- frag_counts %>%
  group_by(library) %>%
  summarise(total = sum(count))
print(counts_by_lib)

# Use non-acrocentric arm positions to count fragments per arm
cat("Counting fragments per chromosome arm...\n")
arm_counts <-
  frag_counts %>%
  # Annotate arms and filter to non-acro arms
  filter(chr %in% autosomes) %>%
  left_join(noacro, by = "chr") %>%
  mutate(arm = ifelse(end < armstart_q, "p", "q")) %>%
  filter(chr %in% noacro_chroms | (chr %in% acro_chroms & arm == "q")) %>%
  group_by(library, chr, arm) %>%
  summarise(count = sum(count)) %>%
  left_join(counts_by_lib, by = "library") %>%
  mutate(prop = count / total)

# Find mean and standard deviation for healthy libraries
cat("Finding mean and standard deviation for healthy libraries...\n")
healthy_arms <-
  arm_counts %>%
  filter(library %in% healthy_libs) %>%
  group_by(chr, arm) %>%
  summarize(mean = mean(prop, na.rm = TRUE), sd = sd(prop, na.rm = TRUE))

# Determine per arm z-scores
cat("Calculating fragment proportion z scores...\n")
arm_z <-
  arm_counts %>%
  left_join(healthy_arms, by = c("chr", "arm")) %>%
  mutate(z = (prop - mean) / sd) %>%
  select(!c("mean", "sd"))

write_tsv(arm_z, file = armz_tsv)
cat("Script completed successfully.\n")