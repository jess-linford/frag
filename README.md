# Fragmentomics Analysis for ctDNA

This repository contains several pipelines and scripts to perform fragmentomic analysis on cfDNA. The frag.smk pipeline will perform complete analysis with a single cutpoint separating short and long fragment lengths. The frag_mult_cutpoints pipeline will perform complete analysis and will allow you to specify a range of cutpoints to separate short and long fragment lengths. It will output separate files for each cutpoint. These two pipelines perform the following analyses:
_ Separate a genome fasta file into windows of a specified width and exclude any blacklist regions.
_ Filter bam files to exclude unmapped or duplicate reads, non-primary alignments, reads with map quality <= 20, and reads longer than a maximum length specified in the pipeline. Optionally, it will filter to reads on autosomes only.
_ Convert bam files to bed files and annotate with GC content.
_ Create fragment length distributions for each library (number of reads of each size) and output in a format that can easily be used for PCA or NMF or graphed by library.
_ Find the median fragment length in each window for each library.
_ Determine the number of fragments in each window and the short to long fragment ratio for each window. Fragment counts are normalized based on the library size and GC content of healthy libraries following the process outlined in Mathios, et al. 2021. https://doi.org/10.1038/s41467-021-24994-w
_ Compute z-scores for the number of fragments on each chromosomal arm compared to healthy libraries as described in Mathios, et al., with the modification that z-score calculations are based on the proportion of fragments on each arm rather than the raw count of fragments on each arm.

Separate pipelines are provided to perform parts of this process in isolation:
_ filter_bam.smk filters bam files as described above or by size only depending on the script specified in the snakemake rule.
_ length_distros_from_bams.smk and length_distros_from_beds.smk will create fragment length distributions from bam or bed files, respectively.
_ make_keep_bed.smk will separate a genome into windows of a specified length and remove blacklist regions.
_ med_window_length.smk will calculate median fragment lengths by window for each library. It requires the bed file that is output from make_keep_bed.smk as an input.
_ ratios_and_z_scores.smk and ratios_and_z_scores_mult_cutpoints.smk will output normalized fragment counts, short and long fragment counts, short:long ratios, and chromosome arm z-scores for either a single short:long cutpoint or a range of cutpoints, respectively.

Parameters:
_ threads = number of available threads to use for a single rule.
_ window_size = width of windows to split the genome into (eg. 5Mb = 5000000).
_ max_filter_length = maximum fragment length to keep in bam and bed files.
_ frag_length_low & frag_length_high = fragment lengths considered in the following calculations: median fragment window lengths, fragment counts, short:long ratios, and chromosome arm z-scores. The fragment length distributions output two sets of files - one showing all fragments below max_filter_length, and one restricted to fragments between frag_length_low and frag_length_high.
_ cutpoint or cutpoints = a single fragment lengths or list of fragment lengths (eg. \[90, 100, 110]) specifying what length separates short and long fragments. Short fragments have lengths between frag_length_low and cutpoint and long fragments have lengths between cutpoint and frag_length_high. The multiple cutpoints pipelines allow you to test a variety of fragment length definitions without running the pipeline more than once.

Input files:
genome_fasta = Fasta reference file for the genome.
blklist = Blacklist regions for the genome. A blacklist for hg38 is available in the ref folder of this repo.
libraries_file = A tab-separated file with at least the following columns: library (name), file (bam file path), and cohort (one cohort should be 'healthy'). Other columns may be included if additional filtering is desired. An example library file is in the ref folder.
cytobands = File with chromosome arm definitions. File for hg38 is in ref folder.

The pipelines are set up to use folders in the following structure, but file paths may be modified to preference.
    parent_directory
    ├── /analysis
    │   ├── /bams
    │   ├── /beds
    │   ├── /length_distros
    │   ├── /medians
    │   ├── /counts
    │   └── /gc_distros
    ├── /benchmark
    ├── /logs
    ├── /ref
    └── /scripts


