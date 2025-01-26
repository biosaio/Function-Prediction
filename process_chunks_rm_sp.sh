#!/bin/bash

mkdir -p blast_output_2

for chunk in blast_output/*.tsv; do
    echo "Processing: $chunk"
    output="blast_output_2/$(basename "$chunk" .tsv).tsv"
    cat "$chunk" | sed 's/sp|//g' | sed 's/|//g' > "$output"
done