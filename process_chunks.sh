#!/bin/bash

mkdir -p blast_output

for chunk in fasta/train_chunks/*.fasta; do
    echo "$chunk"
    output="blast_output/$(basename $chunk .fasta).tsv"
    blastp -db blastdb -query "$chunk" -outfmt "6 qseqid sseqid bitscore evalue score" -out "$output" -num_threads 4
done