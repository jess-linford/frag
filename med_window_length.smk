import pandas as pd
import numpy as np

# Parameters
threads = 5
# Fragment length range
frag_length_low = 80
frag_length_high = 220

# Directory values
parentdir               = "/aclm350-zpool1/jlinford/frag"
analysis_dir            = parentdir + "/analysis"
beds_dir                = parentdir + "/analysis/beds"
medians_dir             = parentdir + "/analysis/medians"
benchdir                = parentdir + "/benchmark"
logdir                  = parentdir + "/logs"
refdir                  = parentdir + "/ref"
scriptdir               = parentdir + "/scripts"

# Input files
keep_bed = refdir + "/keep_5mb.bed"

rule all:
    input:
        expand(medians_dir + "/{library}_med_frag_window_lengths.tsv", library = wildcards.library),
        analysis_dir + "/med_frag_window_lengths.tsv"

# Calculate median window lengths for each library
rule med_frag_window_lengths:
    benchmark: benchdir + "/{library}_med_frag_window_lengths.benchmark.txt",
    input: 
        frag_bed = beds_dir + "/{library}_frag.bed",
        keep_bed = keep_bed,
    log: logdir + "/{library}_med_frag_window_lengths.log",
    output: medians_dir + "/{library}_med_frag_window_lengths.tsv",
    params:
        script = scriptdir + "/med_frag_window_lengths.R",
        frag_length_low = frag_length_low,
        frag_length_high = frag_length_high,
        threads = threads,
    shell:
        """
        Rscript {params.script} \
        {input.frag_bed} \
        {input.keep_bed} \
        {params.frag_length_low} \
        {params.frag_length_high} \
        {params.threads} \
        {output} \
        {log} \
        > {log} 2>&1
        """

# Combine median window lengths for all libraries into one file
rule median_merge:
    benchmark: benchdir + "/median_merge.benchmark.txt",
    input: expand(medians_dir + "/{library}_med_frag_window_lengths.tsv", library = wildcards.library),
    log: logdir + "/median_merge.log",
    output: analysis_dir + "/med_frag_window_lengths.tsv",
    params:
        medians_dir = medians_dir,
        script = scriptdir + "/median_merge.sh",
    shell:
        """
        {params.script} \
        {params.medians_dir} \
        {output} &> {log}
        """