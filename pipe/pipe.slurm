#!/usr/bin/env bash
#SCHEDULER="$HOME/fmacrae/megadaph.private/megadaph/pipe/sbatch_scheduler/sbatch_scheduler.py"
# snakemake --keep-going --cluster "$SCHEDULER" -j 1000 --rerun-incomplete \
#     --restart-times 5 --latency-wait 20 "$@" > snakemake.log 2>&1

snakemake --keep-going -j 1000 --rerun-incomplete \
    --restart-times 5 --latency-wait 20 "$@" > snakemake.log 2>&1
