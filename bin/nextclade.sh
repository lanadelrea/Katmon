#!/usr/bin/env bash

# installs latest nextclade build and says yes to prompts

fasta_path=$1
output_dir=$2

nextclade dataset get --name 'sars-cov-2' --output-dir "${output_dir}/data/sars-cov-2"

nextclade run \
    --input-dataset "${output_dir}/data/sars-cov-2" \
    --output-tsv=nextclade.tsv \
    ${fasta_path}
