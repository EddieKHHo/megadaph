#!/usr/bin/env bash
snakemake --cluster 'sbatch -A schaack --mem=10G --tasks-per-node {threads}' \
    -j "$@"
