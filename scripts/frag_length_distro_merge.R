#!/usr/bin/env Rscript

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
# args[1] is a space-separated list of fragment count file names. 
# frag_length_distro_files is a vector of file names.
frag_length_distro_files <- unlist(strsplit(args[1], " ")) 
output_file_long <- args[2] #output
output_file_wide <- args[3] #output
output_file_long_filtered <- args[4] #output
output_file_wide_filtered <- args[5] #output
frag_length_low <- as.integer(args[6]) #parameter
frag_length_high <- as.integer(args[7]) #parameter
num_threads <- as.integer(args[8]) #parameter
log_file <- args[9]

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

# Function to read and process each file
read_and_process <- function(file) {
  df <- fread(file, col.names = c("length", "count"))
  df[, library := sub("_frag_size_distro\\.tsv$", "", basename(file))]  # Extract the part before "_frag_size_distro.tsv"
  return(df)
}

# Read and concatenate all files using parallel processing
cat("Reading and processing fragment length distribution files...\n")
data_list <- lapply(frag_length_distro_files, read_and_process)
concatenated_data <- rbindlist(data_list)

# Ensure the columns are in the correct order
concatenated_data <- concatenated_data[, .(library, length, count)]

# Write the concatenated data to the long format file
cat("Writing concatenated data to file...\n")
fwrite(concatenated_data, output_file_long)

# Filter the data based on frag_length_low and frag_length_high
cat("Filtering data to desired fragment length range...\n")
filtered_data <- concatenated_data[length >= frag_length_low & length <= frag_length_high]

# Write the filtered concatenated data to the long format file
cat("Writing filtered data to file...\n")
fwrite(filtered_data, output_file_long_filtered)

# Transform data to wide format
cat("Transforming data to wide format...\n")
wide_data <- dcast(concatenated_data, library ~ length, value.var = "count", fill = 0)

# Write the wide format data to file
cat("Writing wide format data to file...\n")
fwrite(wide_data, output_file_wide)

# Transform filtered data to wide format
cat("Transforming filtered data to wide format...\n")
wide_data_filtered <- dcast(filtered_data, library ~ length, value.var = "count", fill = 0)

# Write the filtered wide format data to file
cat("Writing filtered wide format data to file...\n")
fwrite(wide_data_filtered, output_file_wide_filtered)

cat("Script completed successfully\n")