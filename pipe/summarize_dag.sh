#!/usr/bin/env bash
snakemake --dryrun --summary | grep pending
