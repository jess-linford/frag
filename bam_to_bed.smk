import pandas as pd
import numpy as np

# Parameters
threads = 5

# Directory values
# Modify parentdir, bams_dir, and beds_dir as needed
parentdir               = "/aclm350-zpool1/jlinford/test/frag_test"
benchdir                = parentdir + "/benchmark"
logdir                  = parentdir + "/logs"
refdir                  = parentdir + "/ref"
scriptdir               = "scripts"
bams_dir                = parentdir + "/analysis/bams"
beds_dir                = parentdir + "/analysis/beds"

# Input files
genome_fasta = refdir + "/GRCh38.p13.genome.fa"
libraries_file = refdir + "/libraries.txt"

# Read libraries column and extract library column as list
libraries = pd.read_table(libraries_file)
ALL_LIBRARIES = libraries['library'].tolist()

rule all:
    input:
        expand(beds_dir + "/{library}_frag.bed", library = ALL_LIBRARIES)

# Make bedfiles from filtered bams
# Make sure naming pattern in input rule matches the naming pattern of the input files
rule bam_to_bed:
    benchmark: benchdir + "/{library}_bam_to_bed.benchmark.txt",
    input: bams_dir + "/{library}_filt.bam",
    log: logdir + "/{library}_bam_to_bed.log",
    output: beds_dir + "/{library}_frag.bed",
    params:
        fasta = genome_fasta,
        script = scriptdir + "/bam_to_bed.sh",
    threads: threads,
    shell:
        """
        {params.script} \
	    {input} \
        {params.fasta} \
        {threads} \
        {output} > {log} 2>&1
        """
        