# Assign input arguments to variables
fasta_file=$1
window_size=$2
output=$3

# Generate the index file using samtools faidx
samtools faidx "$fasta_file"

# Run the commands with the provided parameters
bedtools makewindows -g "${fasta_file}.fai" -w "$window_size" | \
bedtools nuc -fi "$fasta_file" -bed stdin | \
awk 'NR>1 {print $1 "\t" $2 "\t" $3 "\t" $5}' > "$output"