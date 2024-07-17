import pandas as pd
import numpy as np

# Parameters
window_size = 5000000

# Directory values
parentdir               = "/aclm350-zpool1/jlinford/test/frag_test"
benchdir                = parentdir + "/benchmark"
logdir                  = parentdir + "/logs"
refdir                  = parentdir + "/ref"
scriptdir               = parentdir + "/scripts"

# Input files
genome_fasta = refdir + "/GRCh38.p13.genome.fa"
blklist = refdir + "/hg38-blacklist.v2.bed.gz"

rule all:
    input:
        refdir + "/gc5mb.bed",
        refdir + "/keep_5mb.bed"

# Make 5mb window bed file with GC content from fasta file
rule make_gc_window_bed:
    benchmark: benchdir + "/make_gc_window_bed.benchmark.txt",
    input:
        fasta = genome_fasta
    log: logdir + "/make_gc_window_bed.log",
    output: refdir + "/gc5mb.bed",
    params:
        script = scriptdir + "/make_gc_window_bed.sh",
        window_size = window_size,
    shell:
        """
        {params.script} \
        {input.fasta} \
        {params.window_size} \
        {output} &> {log}
        """

# Make GC and mappability restricted bins
rule make_keep_bed:
    benchmark: benchdir + "/make_keep_bed.benchmark.txt",
    input:
        gc_window_bed = refdir + "/gc5mb.bed",
        blklist = blklist,
    log: logdir + "/make_keep_bed.log",
    output: refdir + "/keep_5mb.bed",
    params:
        script = scriptdir + "/make_keep_bed.sh",
    shell:
        """
        {params.script} \
        {input.gc_window_bed} \
        {input.blklist} \
        {output} &> {log}
        """
