#!/usr/bin/env bash
SCHEDULER="$HOME/fmacrae/megadaph.private/megadaph/pipe/sbatch_scheduler.py"
snakemake --cluster "$SCHEDULER" --cluster-config "sbatch_config.yml" -j 499 \
    --rerun-incomplete "$@"
