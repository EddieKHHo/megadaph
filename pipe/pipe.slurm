#!/usr/bin/env bash
SCHEDULER="$HOME/fmacrae/megadaph.private/megadaph/pipe/sbatch_scheduler/sbatch_scheduler.py"
snakemake --cluster "$SCHEDULER" -j 499 --keep-going --immediate-submit \
    --notemp --rerun-incomplete "$@"
