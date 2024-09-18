import pandas as pd
import numpy as np

# Parameters
threads = 5

# Directory values
input_dir               = "/acl45d1/ris/Active/work/pradeep_project/After_alignment/Arpit_WGS_batch2"
output_dir              = "/aclm350-zpool1/jlinford/frag/analysis/bams"
script_dir              = "/aclm350-zpool1/jlinford/frag/scripts"
log_dir                 = "/aclm350-zpool1/jlinford/frag/logs"

# Determine library names from file names
LIBS, = glob_wildcards(input_dir + "/dedup_Aligned_{lib}.bam")

rule all:
    input:
        expand(output_dir + "/{lib}_filt.bam", lib = LIBS)

# Filter bams
# Choose which script to use based on what filtering you want to do
rule filter_bams:
    input: input_dir + "/dedup_Aligned_{lib}.bam",
    output: output_dir + "/{lib}_filt.bam",
    log: log_dir + "/{lib}_filter_bams.log",
    params:
        # Filter out duplicates, unmapped reads, low-quality reads
        script = script_dir + "/filter_bam_no_size.sh",
        threads = threads,
    shell:
        """
        {params.script} \
        {input} \
        {params.threads} \
        {output} &> {log}
        """
        