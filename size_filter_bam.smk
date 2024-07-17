import pandas as pd
import numpy as np

# Parameters
threads = 5

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
        script = scriptdir + "/size_filter_bam.sh",
        threads = threads,
    shell:
        """
        {params.script} \
        {input} \
        {params.threads} \
        {output} &> {log}
        """