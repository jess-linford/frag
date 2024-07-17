import pandas as pd
import numpy as np

# Parameters
threads = 5
max_filter_length = 1000

# Directory values
parentdir               = "/aclm350-zpool1/jlinford/test/frag_test"
bams_dir                = parentdir + "/analysis/bams"
benchdir                = parentdir + "/benchmark"
logdir                  = parentdir + "/logs"
refdir                  = parentdir + "/ref"
scriptdir               = parentdir + "/scripts"

libraries_file = refdir + "/libraries.tsv"

# Read libraries column and extract library column as list
libraries = pd.read_table(libraries_file)
ALL_LIBRARIES = libraries['library'].tolist()

rule all:
    input:
        expand(bams_dir + "/{library}_filt.bam", library = ALL_LIBRARIES)

# Filter bams
# Choose which script to use based on what filtering you want to do
rule filter_bams:
    benchmark: benchdir + "/{library}_filter_bams.benchmark.txt",
    input: bams_dir + "/{library}.bam",
    log: logdir + "/{library}_filter_bams.log",
    output: bams_dir + "/{library}_filt.bam",
    params:
        # Filter out duplicates, unmapped reads, low-quality reads, and reads <= max_filter_length
        script = scriptdir + "/filter_bam.sh",
        # # Filter to only autosomes
        # script = scriptdir + "/filter_bam_autosomes.sh", 
        # # Filter by size only
        # script = scriptdir + "/size_filter_bam.sh",
        threads = threads,
        max_filter_length = max_filter_length,
    shell:
        """
        {params.script} \
        {input} \
        {params.threads} \
        {params.max_filter_length} \
        {output} &> {log}
        """
        