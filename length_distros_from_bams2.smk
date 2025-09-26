import pandas as pd
import numpy as np

# Parameters
threads = 8
max_filter_length = 1000 # Maximum fragment length to keep in histograms
# Fragment lengths to use in filtered distro files
frag_length_low = 80
frag_length_high = 220

# Directory values
# Modify parentdir, bams_dir, and length_distros_dir as needed
parentdir               = "/aclm350-zpool1/jlinford/frag"
analysis_dir            = parentdir + "/analysis"
bams_dir                = parentdir + "/analysis/bams"
length_distros_dir      = parentdir + "/analysis/length_distros"
benchdir                = parentdir + "/benchmark"
logdir                  = parentdir + "/logs"
scriptdir               = parentdir + "/scripts"

# Determine library names from file names
ALL_LIBRARIES, = glob_wildcards(bams_dir + "/{library}_filt.bam")

rule all:
    input:
        expand(length_distros_dir + "/{library}_frag_length_distro.tsv", library = ALL_LIBRARIES),
        expand(length_distros_dir + "/{library}_frag_length_histogram.pdf", library = ALL_LIBRARIES),
        analysis_dir + "/frag_length_distros_long.tsv",
        analysis_dir + "/frag_length_distros_wide.tsv",
        analysis_dir + "/frag_length_distros_long_filtered.tsv",
        analysis_dir + "/frag_length_distros_wide_filtered.tsv"

# Generate fragment length distributions
rule frag_length_distro:
    benchmark: benchdir + "/{library}_frag_length_distro.benchmark.txt",
    input: bams_dir + "/{library}_filt.bam",
    log: logdir + "/{library}_frag_length_distro.log",
    output: 
        metrics = length_distros_dir + "/{library}_frag_length_distro.tsv",
        histogram = length_distros_dir + "/{library}_frag_length_histogram.pdf",
    params:
        script = scriptdir + "/frag_length_distro_from_bam.sh",
        max_length = max_filter_length,
        outdir = length_distros_dir,
    shell:
        """
        {params.script} \
        {params.outdir} \
        {input} \
        {output.metrics} \
        {output.histogram} \
        {params.max_length} > {log} 2>&1
        """

# Merge fragment length distribution files
rule frag_length_distro_merge:
    benchmark: benchdir + "/frag_length_distro_merge.benchmark.txt",
    input: expand(length_distros_dir + "/{library}_frag_length_distro.tsv", library = ALL_LIBRARIES),
    log: logdir + "/frag_length_distro_merge.log",
    output: 
        long = analysis_dir + "/frag_length_distros_long.tsv",
        wide = analysis_dir + "/frag_length_distros_wide.tsv",
    params:
        script = scriptdir + "/frag_length_distro_merge.R",
    threads: threads,
    shell:
        """
        Rscript {params.script} \
        "{input}" \
        {output.long} \
        {output.wide} \
        {threads} \
        {log} > {log} 2>&1
        """
