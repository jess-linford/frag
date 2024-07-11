start=$SECONDS

snakemake -s frag.smk -j 30 -w 60 -p

end=$SECONDS
duration=$(( end - start ))
echo "frag pipeline took $duration seconds to complete"