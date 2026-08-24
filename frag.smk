#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#     Fragmentomic Analysis of Cell-free DNA Whole Genome Sequencing           #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8
import os
import pandas as pd
import numpy as np

# Load configuration
configfile: "config/config_frag.yaml"

# Parameters from config
threads = config["threads"]
window_size = config["window_size"]
max_filter_length = config["max_filter_length"]
frag_length_low = config["frag_length"]["low"]
frag_length_high = config["frag_length"]["high"]
cutpoint = config["frag_length"]["cutpoint"]
bam_name_pattern = config["bam_name_pattern"]
motif_length = config["motif_length"]

# Directory roots from config
OUTPUT_DIR = config["directories"]["output_folder"]
input_bams_dir = config["directories"]["input_bams"]
refdir = config["directories"]["reference"]

# Internal output directory structure
# These are workflow implementation details rather than user-configurable paths.
DIRS = {
    "filtered_bams": os.path.join(OUTPUT_DIR, "filtered_bams"),
    "beds": os.path.join(OUTPUT_DIR, "beds"),
    "length_distros": os.path.join(OUTPUT_DIR, "length_distros"),
    "medians": os.path.join(OUTPUT_DIR, "medians"),
    "counts": os.path.join(OUTPUT_DIR, "counts"),
    "gc_distros": os.path.join(OUTPUT_DIR, "gc_distros"),
    "end_motifs": os.path.join(OUTPUT_DIR, "end_motifs"),
    "reference_data": os.path.join(OUTPUT_DIR, "reference_data"),
    "bench": os.path.join(OUTPUT_DIR, "benchmark"),
    "logs": os.path.join(OUTPUT_DIR, "logs"),
}

filtered_bams_dir = DIRS["filtered_bams"]
beds_dir = DIRS["beds"]
length_distros_dir = DIRS["length_distros"]
medians_dir = DIRS["medians"]
counts_dir = DIRS["counts"]
gc_distros_dir = DIRS["gc_distros"]
end_motifs_dir = DIRS["end_motifs"]
derived_ref_dir = DIRS["reference_data"]
benchdir = DIRS["bench"]
logdir = DIRS["logs"]

# Derived reference resources are pipeline outputs, not inputs.
gc_window_bed = os.path.join(derived_ref_dir, "gc5mb.bed")
keep_bed = os.path.join(derived_ref_dir, "keep_5mb.bed")

# Scripts live in scripts/ alongside this workflow's Snakefile.
scriptdir = os.path.join(workflow.basedir, "scripts")

# Reference filenames are relative to the configured reference directory.
genome_fasta = os.path.join(refdir, config["references"]["genome_fasta"])
blklist = os.path.join(refdir, config["references"]["blacklist"])
libraries_file = os.path.join(refdir, config["references"]["libraries_file"])
cytobands = os.path.join(refdir, config["references"]["cytobands"])

# Construct input BAM paths from the configured input directory and naming pattern.
# This keeps sample discovery/readability checks consistent with rule filter_bams.
def get_input_bam(library):
    return os.path.join(
        input_bams_dir,
        bam_name_pattern.format(library=library)
    )

# Setup sample name index as a pandas DataFrame.
libraries = pd.read_table(libraries_file)

# Keep only libraries whose configured input BAM is readable.
libraries["input_bam"] = libraries["library"].apply(get_input_bam)
libraries["readable"] = libraries["input_bam"].apply(lambda path: os.access(path, os.R_OK))
libraries = libraries[libraries["readable"]]

ALL_LIBRARIES = libraries['library'].tolist()

# Make a list of healthy libraries
HEALTHY_LIBRARIES = libraries[libraries['cohort'] == 'healthy']['library'].tolist()

# Ensure the output tree exists before any jobs start. This includes directories
# used only for logs and benchmarks as well as scientific outputs.
onstart:
    for directory in [OUTPUT_DIR, *DIRS.values()]:
        os.makedirs(directory, exist_ok=True)

rule all:
    input:
        gc_window_bed,
        keep_bed,
        expand(filtered_bams_dir + "/{library}_filt.bam", library = ALL_LIBRARIES),
        expand(beds_dir + "/{library}_frag.bed", library = ALL_LIBRARIES),
        expand(length_distros_dir + "/{library}_frag_length_distro.tsv", library = ALL_LIBRARIES),
        OUTPUT_DIR + "/frag_length_distros_long.tsv",
        OUTPUT_DIR + "/frag_length_distros_wide.tsv",
        expand(medians_dir + "/{library}_med_frag_window_lengths.tsv", library = ALL_LIBRARIES),
        OUTPUT_DIR + "/med_frag_window_lengths_long.tsv",
        OUTPUT_DIR + "/med_frag_window_lengths_wide.tsv",
        expand(end_motifs_dir + "/{lib}_motifs.tsv", lib = ALL_LIBRARIES),
        OUTPUT_DIR + "/motif_counts.tsv",
        OUTPUT_DIR + "/motifs_rel_freq_wide.tsv",
        expand(gc_distros_dir + "/{library}_gc_distro.csv", library = HEALTHY_LIBRARIES),
        gc_distros_dir + "/healthy_med_gc_distro.rds",
        expand(beds_dir + "/{library}_weights.bed", library = ALL_LIBRARIES),
        expand(beds_dir + "/{library}_short_weights.bed", library = ALL_LIBRARIES),
        expand(beds_dir + "/{library}_long_weights.bed", library = ALL_LIBRARIES),
        expand(counts_dir + "/{library}_count_total.tsv", library = ALL_LIBRARIES),
        expand(counts_dir + "/{library}_count_short.tsv", library = ALL_LIBRARIES),
        expand(counts_dir + "/{library}_count_long.tsv", library = ALL_LIBRARIES),
        OUTPUT_DIR + "/frag_counts.tsv",
        OUTPUT_DIR + "/frag_counts_by_len_class.tsv",
        OUTPUT_DIR + "/ratios_long.tsv",
        OUTPUT_DIR + "/ratios_wide.tsv",
        OUTPUT_DIR + "/armz_long.tsv",
        OUTPUT_DIR + "/armz_wide.tsv"

# Make 5mb window bed file with GC content from fasta file
rule make_gc_window_bed:
    benchmark: benchdir + "/make_gc_window_bed.benchmark.txt",
    input:
        fasta = genome_fasta,
    log: logdir + "/make_gc_window_bed.log",
    output: gc_window_bed,
    params:
        script = scriptdir + "/make_gc_window_bed.sh",
        window_size = window_size,
    shell:
        """
        {params.script} \
        {input.fasta} \
        {params.window_size} \
        {output} > {log} 2>&1
        """

# Make GC and mappability restricted bins
rule make_keep_bed:
    benchmark: benchdir + "/make_keep_bed.benchmark.txt",
    input:
        gc_window_bed = gc_window_bed,
        blklist = blklist,
    log: logdir + "/make_keep_bed.log",
    output: keep_bed,
    params:
        script = scriptdir + "/make_keep_bed.sh",
    shell:
        """
        {params.script} \
        {input.gc_window_bed} \
        {input.blklist} \
        {output} > {log} 2>&1
        """

# Filter bams
rule filter_bams:
    benchmark: benchdir + "/{library}_filter_bams.benchmark.txt",
    input:
        lambda wildcards: get_input_bam(wildcards.library),
    log: logdir + "/{library}_filter_bams.log",
    output: filtered_bams_dir + "/{library}_filt.bam",
    params:
        script = scriptdir + "/filter_bam.sh",
        # script = scriptdir + "/filter_bam_autosomes.sh",      # To filter to only autosomes
        max_filter_length = max_filter_length,
    threads: threads,
    shell:
        """
        {params.script} \
        {input} \
        {threads} \
        {params.max_filter_length} \
        {output} > {log} 2>&1
        """

# Make bedfiles from filtered bams
rule bam_to_bed:
    benchmark: benchdir + "/{library}_bam_to_bed.benchmark.txt",
    input: filtered_bams_dir + "/{library}_filt.bam",
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
        long = OUTPUT_DIR + "/frag_length_distros_long.tsv",
        wide = OUTPUT_DIR + "/frag_length_distros_wide.tsv",
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

# Calculate median window lengths for each library
rule med_frag_window_lengths:
    benchmark: benchdir + "/{library}_med_frag_window_lengths.benchmark.txt",
    input: 
        frag_bed = beds_dir + "/{library}_frag.bed",
        keep_bed = keep_bed,
    log: logdir + "/{library}_med_frag_window_lengths.log",
    output: medians_dir + "/{library}_med_frag_window_lengths.tsv",
    params:
        script = scriptdir + "/med_frag_window_lengths.R",
        frag_length_low = frag_length_low,
        frag_length_high = frag_length_high,
    threads: threads,
    shell:
        """
        Rscript {params.script} \
        {input.frag_bed} \
        {input.keep_bed} \
        {params.frag_length_low} \
        {params.frag_length_high} \
        {threads} \
        {output} \
        {log} > {log} 2>&1
        """

# Combine median window lengths for all libraries into one file
rule median_merge:
    benchmark: benchdir + "/median_merge.benchmark.txt",
    input: expand(medians_dir + "/{library}_med_frag_window_lengths.tsv", library = ALL_LIBRARIES),
    log: logdir + "/median_merge.log",
    output: 
        medians_long = OUTPUT_DIR + "/med_frag_window_lengths_long.tsv",
        medians_wide = OUTPUT_DIR + "/med_frag_window_lengths_wide.tsv",
    params:
        medians_dir = medians_dir,
        script = scriptdir + "/median_merge.sh",
    shell:
        """
        {params.script} \
        {params.medians_dir} \
        {output.medians_long} \
        {output.medians_wide} > {log} 2>&1
        """

# Create end motif files from filtered bams
rule process_motifs:
    benchmark: benchdir + "/{lib}_process_motifs.benchmark.txt",
    input: filtered_bams_dir + "/{lib}_filt.bam",
    log: logdir + "/{lib}_process_motifs.log",
    output: end_motifs_dir + "/{lib}_motifs.tsv",
    params:
        script = scriptdir + "/process_motifs.sh",
        fasta = genome_fasta,
        motif_length = motif_length,
    threads: threads,
    shell:
        """
        {params.script} \
        {input} \
        {params.fasta} \
        {params.motif_length} \
        {threads} \
        {output} > {log} 2>&1
        """

# Aggregate motif files into count matrix and relative frequency matrix
rule count_motifs:
    benchmark: benchdir + "/count_motifs.benchmark.txt",
    input: expand(end_motifs_dir + "/{lib}_motifs.tsv", lib = ALL_LIBRARIES),
    log: logdir + "/count_motifs.log",
    output:
        counts  = OUTPUT_DIR + "/motif_counts.tsv",
        relfreq = OUTPUT_DIR + "/motifs_rel_freq_wide.tsv",
    params:
        script = scriptdir + "/motif_frequency_counter_efficient.py",
        input_dir = end_motifs_dir,
        motif_length = motif_length,
        output_dir = OUTPUT_DIR,
    threads: threads,
    shell:
        """
        python {params.script} \
        {params.input_dir} \
        {params.motif_length} \
        {threads} \
        {params.output_dir} > {log} 2>&1
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
    threads: threads,
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output} \
        {params.frag_length_low} \
        {params.frag_length_high} \
        {threads} \
        {log} > {log} 2>&1
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
        {log} > {log} 2>&1
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
        frag_length_low = frag_length_low,
        cutpoint = cutpoint,
        frag_length_high = frag_length_high,
    threads: threads,
    retries: 5,
    shell:
        """
        Rscript {params.script} \
        {input.healthy_med} \
        {input.frag_bed} \
        {threads} \
        {params.frag_length_low} \
        {params.cutpoint} \
        {params.frag_length_high} \
        {output.weights_bed} \
        {output.short} \
        {output.long} \
        {log} > {log} 2>&1
        """

# Find "counts" (sum of weights) of short and long fragments in each 5mb window
rule frag_window_count:
    benchmark: benchdir + "/{library}_frag_window_count.benchmark.txt",
    input:
        weights_bed = beds_dir + "/{library}_weights.bed",
        short = beds_dir + "/{library}_short_weights.bed",
        long = beds_dir + "/{library}_long_weights.bed",
        keep_bed = keep_bed,
    log: logdir + "/{library}_frag_window_count.log",
    output:
        total_counts = counts_dir + "/{library}_count_total.tsv",
        short = counts_dir + "/{library}_count_short.tsv",
        long = counts_dir + "/{library}_count_long.tsv",
    params:
        script = scriptdir + "/frag_window_count.R",
    threads: threads,
    retries: 5,
    shell:
        """
        Rscript {params.script} \
        {input.weights_bed} \
        {input.keep_bed} \
        {output.total_counts} \
        {threads} \
        {log} > {log} 2>&1
        Rscript {params.script} \
        {input.short} \
        {input.keep_bed} \
        {output.short} \
        {threads} \
        {log} >> {log} 2>&1
        Rscript {params.script} \
        {input.long} \
        {input.keep_bed} \
        {output.long} \
        {threads} \
        {log} >> {log} 2>&1
        """

# Merge fragment counts files from all libraries
# frag_counts.tsv gives total counts by window for each library
# frag_counts_by_len_class.tsv gives short and long counts by window for each library
rule count_merge:
    benchmark: benchdir + "/count_merge.benchmark.txt",
    input: expand(counts_dir + "/{library}_count_{length}.tsv", library = ALL_LIBRARIES, length = ["short", "long", "total"]),
    log: logdir + "/count_merge.log",
    output: 
        frag_counts = OUTPUT_DIR + "/frag_counts.tsv",
        frag_counts_by_len_class = OUTPUT_DIR + "/frag_counts_by_len_class.tsv",
    params:
        script = scriptdir + "/count_merge.R",
    threads: threads,        
    shell:
        """
        Rscript {params.script} \
        "{input}" \
        {output.frag_counts} \
        {output.frag_counts_by_len_class} \
        {threads} \
        {log} > {log} 2>&1
        """

# Calculate short:long fragment ratios and center at 0
rule make_ratios:
    benchmark: benchdir + "/make_ratios.benchmark.txt",
    input: OUTPUT_DIR + "/frag_counts_by_len_class.tsv",
    log: logdir + "/make_ratios.log",
    output: 
        ratios_long = OUTPUT_DIR + "/ratios_long.tsv",
        ratios_wide = OUTPUT_DIR + "/ratios_wide.tsv",
    params:
        script = scriptdir + "/make_ratios.R",
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output.ratios_long} \
        {output.ratios_wide} \
        {log} > {log} 2>&1
        """

# Calculate chromosome arm fragment count z-scores
# If you don't care about duplicate healthy libraries, comment out line 44 of script
rule arm_z_scores:
    benchmark: benchdir + "/arm_z_scores.benchmark.txt",
    input: 
        frag_counts = OUTPUT_DIR + "/frag_counts.tsv",
        cytobands = cytobands,
        libraries = libraries_file,
    log: logdir + "/arm_z_scores.log",
    output: 
        armz_long = OUTPUT_DIR + "/armz_long.tsv",
        armz_wide = OUTPUT_DIR + "/armz_wide.tsv",
    params: 
        script = scriptdir + "/arm_z.R",
    shell:
        """
        Rscript {params.script} \
        {input.cytobands} \
        {input.libraries} \
        {input.frag_counts} \
        {output.armz_long} \
        {output.armz_wide} \
        {log} > {log} 2>&1
        """
