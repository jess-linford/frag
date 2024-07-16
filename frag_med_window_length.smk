#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#     Fragmentomic Analysis of Cell-free DNA Whole Genome Sequencing           #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8
import pandas as pd
import re
import numpy as np
import os

# Parameters
threads = 5
window_size = 5000000
frag_length_low = 30
frag_length_high = 500

# Directory values
datadir                 = "/aclm350-zpool1/jlinford/frag"
frag                    = datadir + "/analysis"
frag_bams               = datadir + "/analysis/bams"
frag_beds               = datadir + "/analysis/beds"
frag_medians            = datadir + "/analysis/medians"
benchdir                = datadir + "/benchmark"
logdir                  = datadir + "/logs"
refdir                  = datadir + "/ref"
scriptdir               = datadir + "/scripts"

# Input files
genome_fasta = refdir + "/GRCh38.p13.genome.fa"
blklist = refdir + "/hg38-blacklist.v2.bed.gz"
libraries_file = refdir + "/libraries.tsv"
cytobands = refdir + "/cytoBand.txt"

# Setup sample name index as a python dictionary
libraries = pd.read_table(libraries_file)

readable = []
for x in libraries.file:
    readable.append(os.access(x, os.R_OK))
# Ensure readable fastqs
libraries['readable'] = readable
libraries = libraries[libraries.readable == True]

# Specify which samples to include -- here excluding serial samples
libraries = libraries[libraries.serial == 0]

# Make the dictionary
library_indict = libraries["library"].tolist()
file_indict = libraries["file"].tolist()
lib_dict = dict(zip(library_indict, file_indict))

ALL_LIBRARIES = list(lib_dict.keys())

rule all:
    input:
        # refdir + "/gc5mb.bed",
        # refdir + "/keep_5mb.bed",
        # expand(frag_bams + "/{library}_filt.bam", library = ALL_LIBRARIES),
        # expand(frag_beds + "/{library}_frag.bed", library = ALL_LIBRARIES),
        expand(frag_medians + "/{library}_med_frag_window_lengths.tsv", library = ALL_LIBRARIES),
        frag + "/med_frag_window_lengths.tsv"

# # Make 5mb window bed file with GC content from fasta file
# rule make_gc_window_bed:
#     benchmark: benchdir + "/make_gc_window_bed.benchmark.txt",
#     input:
#         fasta = genome_fasta
#     log: logdir + "/make_gc_window_bed.log",
#     output: refdir + "/gc5mb.bed",
#     params:
#         script = scriptdir + "/make_gc_window_bed.sh",
#         window_size = window_size,
#     shell:
#         """
#         {params.script} \
#         {input.fasta} \
#         {params.window_size} \
#         {output} &> {log}
#         """

# # Make GC and mappability restricted bins
# rule make_keep_bed:
#     benchmark: benchdir + "/make_keep_bed.benchmark.txt",
#     input:
#         gc_window_bed = refdir + "/gc5mb.bed",
#         blklist = blklist,
#     log: logdir + "/make_keep_bed.log",
#     output: refdir + "/keep_5mb.bed",
#     params:
#         script = scriptdir + "/make_keep_bed.sh",
#     shell:
#         """
#         {params.script} \
#         {input.gc_window_bed} \
#         {input.blklist} \
#         {output} &> {log}
#         """

# # Filter bams
# rule filter_bams:
#     benchmark: benchdir + "/{library}_filter_bams.benchmark.txt",
#     input: frag_bams + "/{library}.bam",
#     log: logdir + "/{library}_filter_bams.log",
#     output: frag_bams + "/{library}_filt.bam",
#     params:
#         script = scriptdir + "/filter_bam.sh",
#         threads = threads,
#     shell:
#         """
#         {params.script} \
#         {input} \
#         {params.threads} \
#         {output} &> {log}
#         """

# # Make bedfiles from filtered bams
# rule filt_bam_to_frag_bed:
#     benchmark: benchdir + "/{library}_filt_bam_to_frag_bed.benchmark.txt",
#     input: frag_bams + "/{library}_filt.bam",
#     log: logdir + "/{library}_filt_bam_to_frag_bed.log",
#     output: frag_beds + "/{library}_frag.bed",
#     params:
#         fasta = genome_fasta,
#         script = scriptdir + "/filt_bam_to_frag_bed.sh",
#         threads = threads,
#     shell:
#         """
#         {params.script} \
# 	    {input} \
#         {params.fasta} \
#         {params.threads} \
#         {output} > {log} 2>&1
#         """

# Calculate median window lengths for each library
rule med_frag_window_lengths:
    benchmark: benchdir + "/{library}_med_frag_window_lengths.benchmark.txt",
    input: 
        frag_bed = frag_beds + "/{library}_frag.bed",
        keep_bed = refdir + "/keep_5mb.bed",
    log: logdir + "/{library}_med_frag_window_lengths.log",
    output: frag_medians + "/{library}_med_frag_window_lengths.tsv",
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
    input: expand(frag_medians + "/{library}_med_frag_window_lengths.tsv", library=ALL_LIBRARIES),
    log: logdir + "/median_merge.log",
    output: frag + "/med_frag_window_lengths.tsv",
    params:
        medians_dir = frag_medians,
        script = scriptdir + "/median_merge.sh",
    shell:
        """
        {params.script} \
        {params.medians_dir} \
        {output} &> {log}
        """