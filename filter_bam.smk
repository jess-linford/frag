import pandas as pd
import numpy as np

# Parameters
threads = 5
max_filter_length = 1000

# Directory values
# Modify parentdir, inputdir, and outputdir as needed
parentdir               = "/aclm350-zpool1/jlinford/frag"
benchdir                = parentdir + "/benchmark"
logdir                  = parentdir + "/logs"
scriptdir               = parentdir + "/scripts"
inputdir                = "/path/to/input/bams"
outputdir               = "/path/to/output/bams"

# Filter bams
# Make sure naming pattern in input rule matches the naming pattern of the input files
rule filter_bams:
    benchmark: benchdir + "/{library}_filter_bams.benchmark.txt",
    input: inputdir + "/{library}.bam",
    log: logdir + "/{library}_filter_bams.log",
    output: outputdir + "/{library}_filt.bam",
    params:
        script = scriptdir + "/filter_bam.sh",
        # script = scriptdir + "/filter_bam_autosomes.sh",      # To filter to only autosomes
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