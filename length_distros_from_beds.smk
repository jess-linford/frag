import pandas as pd
import numpy as np

# Parameters
threads = 8
max_filter_length = 1000 # Maximum fragment length to keep in histograms

# Directory values
# Modify parentdir, beds_dir, and length_distros_dir as needed
parentdir               = "/aclm350-zpool1/jlinford/test/frag_test"
analysis_dir            = parentdir + "/analysis"
beds_dir                = parentdir + "/analysis/beds"
length_distros_dir      = parentdir + "/analysis/length_distros"
benchdir                = parentdir + "/benchmark"
logdir                  = parentdir + "/logs"
refdir                  = parentdir + "/ref"
scriptdir               = parentdir + "/scripts"

# Input files
libraries_file = refdir + "/libraries.txt"

# Read libraries column and extract library column as list
libraries = pd.read_table(libraries_file)
ALL_LIBRARIES = libraries['library'].tolist()

rule all:
    input:
        expand(length_distros_dir + "/{library}_frag_length_distro.tsv", library = ALL_LIBRARIES),
        analysis_dir + "/frag_length_distros_long.tsv",
        analysis_dir + "/frag_length_distros_wide.tsv",
        analysis_dir + "/frag_length_distros_long_filtered.tsv",
        analysis_dir + "/frag_length_distros_wide_filtered.tsv"

# Generate fragment length distributions
rule frag_length_distro:
    benchmark: benchdir + "/{library}_frag_length_distro.benchmark.txt",
    input: beds_dir + "/{library}_frag.bed",
    log: logdir + "/{library}_frag_length_distro.log",
    output: length_distros_dir + "/{library}_frag_length_distro.tsv",
    params:
        script = scriptdir + "/frag_length_distro_from_bed.R",
        max_length = max_filter_length,
    threads: threads,
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output} \
        {params.max_length} \
        {threads} \
        {log} > {log} 2>&1
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
