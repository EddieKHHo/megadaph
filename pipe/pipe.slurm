#!/usr/bin/env bash
SCHEDULER="$HOME/fmacrae/megadaph.private/megadaph/pipe/sbatch_scheduler/sbatch_scheduler.py"
snakemake --keep-going --cluster "$SCHEDULER" -j 499 --rerun-incomplete \
    --latency-wait 20 "$@"
