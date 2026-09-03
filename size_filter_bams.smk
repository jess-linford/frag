import pandas as pd
import numpy as np

# Parameters
n_threads = 8
min_filter_length = 90
max_filter_length = 150

# Directory values
input_dir  = "/mnt/data/prostate/bams"
output_dir = "/mnt/data/prostate/bams_90_150"
script_dir = "/mnt/data/scripts/frag/scripts"

# Determine library names from file names
LIBS, = glob_wildcards(
    input_dir + "/{lib}.GRCh38.bwa_meth.coorsort.filt.bam"
)

rule all:
    input:
        expand(
            output_dir + "/{lib}_filt_90_150.bam",
            lib=LIBS
        )

# Filter BAMs
rule filter_bams:
    input:
        bam=input_dir + "/{lib}.GRCh38.bwa_meth.coorsort.filt.bam"
    output:
        bam=output_dir + "/{lib}_filt_90_150.bam",
        bai=output_dir + "/{lib}_filt_90_150.bam.bai"
    params:
        script=script_dir + "/size_filter_bam2.sh",
        min_filter_length=min_filter_length,
        max_filter_length=max_filter_length
    threads:
        n_threads
    log:
        output_dir + "/logs/{lib}_filter_bam.log"
    shell:
        """
        mkdir -p $(dirname {log})

        {params.script} \
            {input.bam} \
            {threads} \
            {params.min_filter_length} \
            {params.max_filter_length} \
            {output.bam} \
            > {log} 2>&1
        """