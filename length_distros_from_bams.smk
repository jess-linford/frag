import pandas as pd
import numpy as np

# Parameters
threads = 5
max_filter_length = 1000 # Maximum fragment length to keep in histograms
# Fragment lengths to use in filtered distro files
frag_length_low = 80
frag_length_high = 220

# Directory values
# Modify parentdir, bams_dir, and length_distros_dir as needed
parentdir               = "/aclm350-zpool1/jlinford/test/frag_test"
analysis_dir            = parentdir + "/analysis"
bams_dir                = parentdir + "/analysis/bams"
length_distros_dir      = parentdir + "/analysis/length_distros"
benchdir                = parentdir + "/benchmark"
logdir                  = parentdir + "/logs"
refdir                  = parentdir + "/ref"
scriptdir               = parentdir + "/scripts"

# Input files
libraries_file = refdir + "/libraries.tsv"

# Read libraries column and extract library column as list
libraries = pd.read_table(libraries_file)
ALL_LIBRARIES = libraries['library'].tolist()

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
        {params.max_length} \
        &> {log}
        """

# Merge fragment length distribution files
rule frag_length_distro_merge:
    benchmark: benchdir + "/frag_length_distro_merge.benchmark.txt",
    input: expand(length_distros_dir + "/{library}_frag_length_distro.tsv", library = ALL_LIBRARIES),
    log: logdir + "/frag_length_distro_merge.log",
    output: 
        long = analysis_dir + "/frag_length_distros_long.tsv",
        wide = analysis_dir + "/frag_length_distros_wide.tsv",
        long_filtered = analysis_dir + "/frag_length_distros_long_filtered.tsv",
        wide_filtered = analysis_dir + "/frag_length_distros_wide_filtered.tsv",
    params:
        script = scriptdir + "/frag_length_distro_merge.R",
        frag_length_low = frag_length_low,
        frag_length_high = frag_length_high,
        threads = threads,
    shell:
        """
        Rscript {params.script} \
        "{input}" \
        {output.long} \
        {output.wide} \
        {output.long_filtered} \
        {output.wide_filtered} \
        {params.frag_length_low} \
        {params.frag_length_high} \
        {params.threads} \
        {log} &> {log}
        """
