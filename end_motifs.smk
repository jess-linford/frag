# Load configuration
configfile: "config_end_motifs.yaml"

# Specify folders and parameters
filtered_bams_dir = config["filtered_bams_dir"]
end_motifs_dir = config["end_motifs_dir"]
analysis_dir = config["analysis_dir"]
scriptdir = config["scriptdir"]
genome_fasta = config["genome_fasta"]
motif_length = config["motif_length"]
threads = config["threads"]

# Determine library names from file names
LIBS, = glob_wildcards(filtered_bams_dir + "/{lib}_filt.bam")

rule all:
    input:
        expand(end_motifs_dir + "/{lib}_motifs.tsv", lib = LIBS),
        analysis_dir + "/motif_counts.tsv",
        analysis_dir + "/motifs_rel_freq_wide.tsv"

# Create end motif files from filtered bams
rule process_motifs:
    benchmark: benchdir + "/{lib}_process_motifs.benchmark.txt",
    input: filtered_bams_dir + "/{lib}_filt.bam",
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
    input: expand(end_motifs_dir + "/{lib}_motifs.tsv", lib = LIBS),
    output:
        counts  = analysis_dir + "/motif_counts.tsv",
        relfreq = analysis_dir + "/motifs_rel_freq_wide.tsv",
    params:
        script = scriptdir + "/motif_frequency_counter_efficient.py",
        input_dir = end_motifs_dir,
        motif_length = motif_length,
        output_dir = analysis_dir,
    threads: threads,
    shell:
        """
        python {params.script} \
        {params.input_dir} \
        {params.motif_length} \
        {threads} \
        {params.output_dir} > {log} 2>&1 
        """