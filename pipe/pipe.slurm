#!/usr/bin/env bash
snakemake --cluster 'sbatch -A schaack --mem=100G --tasks-per-node {threads}' \
    -j 384 --rerun-incomplete "$@"
