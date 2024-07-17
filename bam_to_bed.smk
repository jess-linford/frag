import pandas as pd
import numpy as np

# Parameters
threads = 5

# Directory values
# Modify parentdir, input_bams_dir, and output_beds_dir as needed
parentdir               = "/aclm350-zpool1/jlinford/frag"
benchdir                = parentdir + "/benchmark"
logdir                  = parentdir + "/logs"
refdir                  = parentdir + "/ref"
scriptdir               = parentdir + "/scripts"
input_bams_dir          = parentdir + "/analysis/bams"
output_beds_dir         = parentdir + "/analysis/beds"

# Input files
genome_fasta = refdir + "/GRCh38.p13.genome.fa"

# Make bedfiles from filtered bams
# Make sure naming pattern in input rule matches the naming pattern of the input files
rule bam_to_bed:
    benchmark: benchdir + "/{library}_bam_to_bed.benchmark.txt",
    input: input_bams_dir + "/{library}_filt.bam",
    log: logdir + "/{library}_bam_to_bed.log",
    output: output_beds_dir + "/{library}_frag.bed",
    params:
        fasta = genome_fasta,
        script = scriptdir + "/bam_to_bed.sh",
        threads = threads,
    shell:
        """
        {params.script} \
	    {input} \
        {params.fasta} \
        {params.threads} \
        {output} > {log} 2>&1
        """