# 5' end-motif workflow
#
# The extraction step uses alignment strand to identify the biological 5' end
# of each read and bedtools getfasta -s to return minus-strand motifs in the
# correct 5' -> 3' orientation.

configfile: "config/config_end_motifs.yaml"


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

filtered_bams_dir = config["filtered_bams_dir"]
end_motifs_dir = config["end_motifs_dir"]
analysis_dir = config["analysis_dir"]
scriptdir = config["scriptdir"]
genome_fasta = config["genome_fasta"]
motif_length = config["motif_length"]
threads = config["threads"]

if not isinstance(motif_length, int) or motif_length < 1:
    raise ValueError("motif_length must be a positive integer.")

# Keep logs and benchmark files with the analysis outputs unless explicitly
# changed here later.
logdir = analysis_dir + "/logs"
benchdir = analysis_dir + "/benchmarks"


# -----------------------------------------------------------------------------
# Samples
# -----------------------------------------------------------------------------

LIBS, = glob_wildcards(filtered_bams_dir + "/{lib}_filt.bam")


rule all:
    input:
        expand(end_motifs_dir + "/{lib}_motifs.tsv", lib=LIBS),
        analysis_dir + "/motif_counts.tsv",
        analysis_dir + "/motifs_rel_freq_wide.tsv",


# -----------------------------------------------------------------------------
# Extract strand-aware 5' motif counts from each filtered BAM.
#
# Output format for each sample:
#   motif<TAB>count
#
# The extractor counts both mates but determines the 5' end from the actual
# alignment strand (+/-), not from R1/R2 identity.
# -----------------------------------------------------------------------------

rule extract_5p_motifs:
    input:
        bam=filtered_bams_dir + "/{lib}_filt.bam",
        fasta=genome_fasta,
    output:
        end_motifs_dir + "/{lib}_motifs.tsv",
    log:
        logdir + "/{lib}_extract_5p_motifs.log",
    benchmark:
        benchdir + "/{lib}_extract_5p_motifs.benchmark.txt",
    params:
        script=scriptdir + "/extract_5p_motifs.sh",
        motif_length=motif_length,
    threads: threads
    shell:
        r"""
        {params.script:q} \
            {input.bam:q} \
            {input.fasta:q} \
            {output:q} \
            {params.motif_length} \
            {threads} \
            > {log:q} 2>&1
        """


# -----------------------------------------------------------------------------
# Combine the per-sample motif count tables into:
#   1. motif_counts.tsv          motifs x libraries
#   2. motifs_rel_freq_wide.tsv libraries x motifs
#
# The Python script understands the new motif<TAB>count extractor output.
# -----------------------------------------------------------------------------

rule count_motifs:
    input:
        expand(end_motifs_dir + "/{lib}_motifs.tsv", lib=LIBS),
    output:
        counts=analysis_dir + "/motif_counts.tsv",
        relfreq=analysis_dir + "/motifs_rel_freq_wide.tsv",
    log:
        logdir + "/count_motifs.log",
    benchmark:
        benchdir + "/count_motifs.benchmark.txt",
    params:
        script=scriptdir + "/motif_frequency_counter_efficient.py",
        input_dir=end_motifs_dir,
        motif_length=motif_length,
        output_dir=analysis_dir,
    threads: threads
    shell:
        r"""
        python {params.script:q} \
            {params.input_dir:q} \
            {params.motif_length} \
            {threads} \
            {params.output_dir:q} \
            > {log:q} 2>&1
        """
