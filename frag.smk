#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#     Fragmentomic Analysis of Cell-free DNA Whole Genome Sequencing           #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8
import pandas as pd
import numpy as np

# Parameters
threads = 5
window_size = 5000000
max_filter_length = 1000 # Maximum fragment length to keep in bam and bed files
# Fragment length definitions used for short and long fragments (low to cutpoint and cutpoint to high)
frag_length_low = 80
frag_length_high = 220
cutpoint = 150

# Directory values
parentdir               = "/aclm350-zpool1/jlinford/test/frag_test"
analysis_dir            = parentdir + "/analysis"
bams_dir                = parentdir + "/analysis/bams"
beds_dir                = parentdir + "/analysis/beds"
length_distros_dir      = parentdir + "/analysis/length_distros"
medians_dir             = parentdir + "/analysis/medians"
counts_dir              = parentdir + "/analysis/counts"
gc_distros_dir          = parentdir + "/analysis/gc_distros"
benchdir                = parentdir + "/benchmark"
logdir                  = parentdir + "/logs"
refdir                  = parentdir + "/ref"
scriptdir               = parentdir + "/scripts"

# Input files
genome_fasta = refdir + "/GRCh38.p13.genome.fa"
blklist = refdir + "/hg38-blacklist.v2.bed.gz"
# Libraries file is a tab-separated file that must include at least the following columns: library, file(bam), and cohort.
# Toggle commenting out or including lines 51, 56, and 59-62 depending on which libraries you want to include.
libraries_file = refdir + "/libraries.tsv"
cytobands = refdir + "/cytoBand.txt"

# Setup sample name index as a python dictionary
libraries = pd.read_table(libraries_file)

# Ensure readable bams
readable = []
for x in libraries.file:
    readable.append(os.access(x, os.R_OK))
libraries['readable'] = readable
libraries = libraries[libraries.readable == True]

# # To exclude serial samples:
# libraries = libraries[libraries.serial == 0]

ALL_LIBRARIES = libraries['library'].tolist()

# Make a list of healthy libraries if you don't care about duplicates
HEALTHY_LIBRARIES = libraries[libraries['cohort'] == 'healthy']['library'].tolist()

# # To exclude duplicate healthy libraries run in multiple batches:
# HEALTHY_LIBRARIES = libraries[
#     (libraries['cohort'] == 'healthy') & 
#     ((libraries['duplicate'] == 0) | ((libraries['duplicate'] == 1) & (libraries['batch'] != 1)))
# ]['library'].tolist()

rule all:
    input:
        refdir + "/gc5mb.bed",
        refdir + "/keep_5mb.bed",
        expand(bams_dir + "/{library}_filt.bam", library = ALL_LIBRARIES),
        expand(beds_dir + "/{library}_frag.bed", library = ALL_LIBRARIES),
        expand(length_distros_dir + "/{library}_frag_length_distro.tsv", library = ALL_LIBRARIES),
        analysis_dir + "/frag_length_distros_long.tsv",
        analysis_dir + "/frag_length_distros_wide.tsv",
        analysis_dir + "/frag_length_distros_long_filtered.tsv",
        analysis_dir + "/frag_length_distros_wide_filtered.tsv",
        expand(medians_dir + "/{library}_med_frag_window_lengths.tsv", library = ALL_LIBRARIES),
        analysis_dir + "/med_frag_window_lengths.tsv",
        expand(gc_distros_dir + "/{library}_gc_distro.csv", library = ALL_LIBRARIES),
        gc_distros_dir + "/healthy_med_gc_distro.rds",
        expand(beds_dir + "/{library}_weights.bed", library = ALL_LIBRARIES),
        expand(beds_dir + "/{library}_short_weights.bed", library = ALL_LIBRARIES),
        expand(beds_dir + "/{library}_long_weights.bed", library = ALL_LIBRARIES),
        expand(counts_dir + "/{library}_count_total.tsv", library = ALL_LIBRARIES),
        expand(counts_dir + "/{library}_count_short.tsv", library = ALL_LIBRARIES),
        expand(counts_dir + "/{library}_count_long.tsv", library = ALL_LIBRARIES),
        analysis_dir + "/frag_counts.tsv",
        analysis_dir + "/frag_counts_by_len_class.tsv",
        analysis_dir + "/ratios.tsv",
        analysis_dir + "/arm_z.tsv"

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

# Filter bams
rule filter_bams:
    benchmark: benchdir + "/{library}_filter_bams.benchmark.txt",
    input: bams_dir + "/{library}.bam",
    log: logdir + "/{library}_filter_bams.log",
    output: bams_dir + "/{library}_filt.bam",
    params:
        script = scriptdir + "/filter_bam.sh",
        # script = scriptdir + "/filter_bam_autosomes.sh",      # To filter to only autosomes
        threads = threads,
        max_filter_length = max_filter_length,
    shell:
        """
        {params.script} \
        {input} \
        {params.threads} \
        {params.max_filter_length} \
        {output} &> {log}
        """

# Make bedfiles from filtered bams
rule bam_to_bed:
    benchmark: benchdir + "/{library}_bam_to_bed.benchmark.txt",
    input: bams_dir + "/{library}_filt.bam",
    log: logdir + "/{library}_bam_to_bed.log",
    output: beds_dir + "/{library}_frag.bed",
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
        {output} &> {log}
        """

# Generate fragment length distributions
rule frag_length_distro:
    benchmark: benchdir + "/{library}_frag_length_distro.benchmark.txt",
    input: beds_dir + "/{library}_frag.bed",
    log: logdir + "/{library}_frag_length_distro.log",
    output: length_distros_dir + "/{library}_frag_length_distro.tsv",
    params:
        script = scriptdir + "/frag_length_distro_from_bed.R",
        max_length = max_filter_length,
        threads = threads,
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output} \
        {params.max_length} \
        {params.threads} \
        {log} &> {log}
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

# Calculate median window lengths for each library
rule med_frag_window_lengths:
    benchmark: benchdir + "/{library}_med_frag_window_lengths.benchmark.txt",
    input: 
        frag_bed = beds_dir + "/{library}_frag.bed",
        keep_bed = refdir + "/keep_5mb.bed",
    log: logdir + "/{library}_med_frag_window_lengths.log",
    output: medians_dir + "/{library}_med_frag_window_lengths.tsv",
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
        {log} &> {log}
        """

# Combine median window lengths for all libraries into one file
rule median_merge:
    benchmark: benchdir + "/median_merge.benchmark.txt",
    input: expand(medians_dir + "/{library}_med_frag_window_lengths.tsv", library = ALL_LIBRARIES),
    log: logdir + "/median_merge.log",
    output: analysis_dir + "/med_frag_window_lengths.tsv",
    params:
        medians_dir = medians_dir,
        script = scriptdir + "/median_merge.sh",
    shell:
        """
        {params.script} \
        {params.medians_dir} \
        {output} &> {log}
        """

# Make GC distributions
rule gc_distro:
    benchmark: benchdir + "/{library}_gc_distro.benchmark.txt",
    input: beds_dir + "/{library}_frag.bed",
    log: logdir + "/{library}_gc_distro.log",
    output: gc_distros_dir + "/{library}_gc_distro.csv",
    params:
        script = scriptdir + "/gc_distro.R",
        frag_length_low = frag_length_low,
        frag_length_high = frag_length_high,
        threads = threads,
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output} \
        {params.frag_length_low} \
        {params.frag_length_high} \
        {params.threads} \
        {log} &> {log}
        """

# Make healthy GC distributions summary file
rule healthy_gc:
    benchmark: benchdir + "/healthy_gc.benchmark.txt",
    input: expand(gc_distros_dir + "/{library}_gc_distro.csv", library = HEALTHY_LIBRARIES),
    log: logdir + "/healthy_gc.log",
    output: gc_distros_dir + "/healthy_med_gc_distro.rds",
    params:
        script = scriptdir + "/healthy_gc.R",
    shell:
        """
        Rscript {params.script} \
        "{input}" \
        {output} \
        {log} &> {log}
        """

# Add weights to fragments to normalize for library size and GC bias and split into short and long fragments
rule frag_weighting_and_size_split:
    benchmark: benchdir + "/{library}_frag_weighting_and_size_split.benchmark.txt",
    input:
        frag_bed = beds_dir + "/{library}_frag.bed",
        healthy_med = gc_distros_dir + "/healthy_med_gc_distro.rds",
    log: logdir + "/{library}_frag_weighting_and_size_split.log",
    output: 
        weights_bed = beds_dir + "/{library}_weights.bed",
        short = beds_dir + "/{library}_short_weights.bed",
        long = beds_dir + "/{library}_long_weights.bed",
    params:
        script = scriptdir + "/frag_weighting_and_size_split.R",
        threads = threads,
        frag_length_low = frag_length_low,
        cutpoint = cutpoint,
        frag_length_high = frag_length_high,
    shell:
        """
        Rscript {params.script} \
        {input.healthy_med} \
        {input.frag_bed} \
        {params.threads} \
        {params.frag_length_low} \
        {params.cutpoint} \
        {params.frag_length_high} \
        {output.weights_bed} \
        {output.short} \
        {output.long} \
        {log} &> {log}
        """

# Find "counts" (sum of weights) of short and long fragments in each 5mb window
rule frag_window_count:
    benchmark: benchdir + "/{library}_frag_window_count.benchmark.txt",
    input:
        weights_bed = beds_dir + "/{library}_weights.bed",
        short = beds_dir + "/{library}_short_weights.bed",
        long = beds_dir + "/{library}_long_weights.bed",
        keep_bed = refdir + "/keep_5mb.bed",
    log: logdir + "/{library}_frag_window_count.log",
    output:
        total_counts = counts_dir + "/{library}_count_total.tsv",
        short = counts_dir + "/{library}_count_short.tsv",
        long = counts_dir + "/{library}_count_long.tsv",
    params:
        script = scriptdir + "/frag_window_count.R",
        threads = threads,
    shell:
        """
        Rscript {params.script} \
        {input.weights_bed} \
        {input.keep_bed} \
        {output.total_counts} \
        {params.threads} \
        {log} &>> {log}
        Rscript {params.script} \
        {input.short} \
        {input.keep_bed} \
        {output.short} \
        {params.threads} \
        {log} &>> {log}
        Rscript {params.script} \
        {input.long} \
        {input.keep_bed} \
        {output.long} \
        {params.threads} \
        {log} &>> {log}
        """

# Merge fragment counts files from all libraries
# frag_counts.tsv gives total counts by window for each library
# frag_counts_by_len_class.tsv gives short and long counts by window for each library
rule count_merge:
    benchmark: benchdir + "/count_merge.benchmark.txt",
    input: expand(counts_dir + "/{library}_count_{length}.tsv", library = ALL_LIBRARIES, length = ["short", "long", "total"]),
    log: logdir + "/count_merge.log",
    output: 
        frag_counts = analysis_dir + "/frag_counts.tsv",
        frag_counts_by_len_class = analysis_dir + "/frag_counts_by_len_class.tsv",
    params:
        script = scriptdir + "/count_merge.R",
        threads = threads,        
    shell:
        """
        Rscript {params.script} \
        "{input}" \
        {output.frag_counts} \
        {output.frag_counts_by_len_class} \
        {params.threads} \
        {log} &> {log}
        """

# Calculate short:long fragment ratios and center at 0
rule make_ratios:
    benchmark: benchdir + "/make_ratios.benchmark.txt",
    input: analysis_dir + "/frag_counts_by_len_class.tsv",
    log: logdir + "/make_ratios.log",
    output: analysis_dir + "/ratios.tsv",
    params:
        script = scriptdir + "/make_ratios.R",
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output} \
        {log} &> {log}
        """

# Calculate chromosome arm fragment count z-scores
# If you don't care about duplicate healthy libraries, comment out line 44 of script
rule arm_z_scores:
    benchmark: benchdir + "/arm_z_scores.benchmark.txt",
    input: 
        frag_counts = analysis_dir + "/frag_counts.tsv",
        cytobands = cytobands,
        libraries = libraries_file,
    log: logdir + "/arm_z_scores.log",
    output: analysis_dir + "/arm_z.tsv",
    params: 
        script = scriptdir + "/arm_z.R",
    shell:
        """
        Rscript {params.script} \
        {input.cytobands} \
        {input.libraries} \
        {input.frag_counts} \
        {output} \
        {log} &> {log}
        """
