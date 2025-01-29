import pandas as pd
import numpy as np

# Parameters
threads = 5
min_filter_length = 90
max_filter_length = 150

# Directory values
input_dir               = "/logo2/irfan/GBM_project/GBM_batch2/FASTQ/aligned_bams/"
output_dir              = "/aclm350-zpool1/jlinford/ichor/GBM_trimmed"
script_dir              = "/aclm350-zpool1/jlinford/frag/scripts"

# Determine library names from file names
LIBS, = glob_wildcards(input_dir + "/{lib}.bam")

rule all:
    input:
        expand(output_dir + "/{lib}_filt.bam", lib = LIBS)

# Filter bams
# Choose which script to use based on what filtering you want to do
rule filter_bams:
    input: input_dir + "/{lib}.bam",
    output: output_dir + "/{lib}_filt.bam",
    params:
        script = script_dir + "/size_filter_bam2.sh",
        threads = threads,
        min_filter_length = min_filter_length,
        max_filter_length = max_filter_length,
    shell:
        """
        {params.script} \
        {input} \
        {params.threads} \
        {params.min_filter_length} \
        {params.max_filter_length} \
        {output}
        """
        