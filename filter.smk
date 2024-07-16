#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#     Fragmentomic Analysis of Cell-free DNA Whole Genome Sequencing           #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8
import pandas as pd
import re
import numpy as np

# Parameters
threads = 5

# Directory values
datadir                 = "/aclm350-zpool1/jlinford/frag"
frag_bams               = datadir + "/analysis/bams"
benchdir                = datadir + "/benchmark"
logdir                  = datadir + "/logs"
refdir                  = datadir + "/ref"
scriptdir               = datadir + "/scripts"

# Input files
libraries_file = refdir + "/libraries.tsv"

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

# Make  a list of healthy libraries -- here excluding duplicates run in multiple batches
HEALTHY_LIBRARIES = libraries[
    (libraries['cohort'] == 'healthy') & 
    ((libraries['duplicate'] == 0) | ((libraries['duplicate'] == 1) & (libraries['batch'] != 1)))
]['library'].tolist()

rule all:
    input:
        expand("/aclm350-zpool1/jlinford/ichor/bams500/{library}_filt_500.bam", library = ALL_LIBRARIES)

# Filter bams
rule filter_bams:
    benchmark: benchdir + "/{library}_filter_bams.benchmark.txt",
    input: frag_bams + "/{library}_filt.bam",
    log: logdir + "/{library}_filter_bams.log",
    output: "/aclm350-zpool1/jlinford/ichor/bams500/{library}_filt_500.bam",
    params:
        script = scriptdir + "/size_filter_bam.sh",
        threads = threads,
    shell:
        """
        echo "Running filter_bams for {input} with {params.threads} threads" >> {log}
        {params.script} \
        {input} \
        {params.threads} \
        {output} &>> {log}
        echo "Finished filter_bams for {output}" >> {log}
        """
