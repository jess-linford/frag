import pandas as pd
import numpy as np

# Parameters
threads = 5

# Directory values
input_dir                = "/home/jupyter/data/noah/analysis/bams"
output_dir               = "/home/jupyter/data/noah/analysis/beds"
script_dir               = "/home/jupyter/data/noah/scripts"

# Input files
genome_fasta = "/home/jupyter/data/noah/ref/GRCh38-DAC-U2AF1.fna"

# Determine library names from file names
LIBRARIES, = glob_wildcards(input_dir + "/{library}_filt.bam")

rule all:
    input:
        expand(output_dir + "/{library}_frag.bed", library = LIBRARIES)

# Make bedfiles from filtered bams
# Make sure naming pattern in input rule matches the naming pattern of the input files
rule bam_to_bed:
    input: input_dir + "/{library}_filt.bam",
    log: output_dir + "/logs/{library}_bam_to_bed.log"
    output: output_dir + "/{library}_frag.bed",
    params:
        fasta = genome_fasta,
        script = script_dir + "/bam_to_bed.sh",
        threads = threads,
    shell:
        """
        {params.script} \
	    {input} \
        {params.fasta} \
        {params.threads} \
        {output} > {log} 2>&1
        """
        