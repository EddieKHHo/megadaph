# Megadaph Pipeline Implementation Notes

The main pipeline script `Snakefile` is implemented in [Snakemake](https://snakemake.readthedocs.io/en/stable/) a self-documenting pipeline development language.

Many external binaries and scripts are required to run the pipeline, most of which are not included in the github repo. The following is an outline of the pipeline implementation for researchers who wish to replicate or expand on these analyses.

## Snakefile

The pipeline is split into four main sections:
1. Read QC
2. *De Novo* Assembly/Decontamination Iterations
3. Assembly reduction/scaffolding
4. Variant Calling

Snakemake turned out to be a poor choice for this pipeline due to its iterative nature. Since Snakemake works backwards from output to input, looped analysis sections can get complicated and messy. The pipeline relies on a lot of buggy getter functions for matching rule outputs -> inputs, which are probably unintelligible (sorry!).

## snakemake_slurm_utils/

This is a small module I wrote to improve Snakemake's job scheduling with [SLURM](https://slurm.schedmd.com/). The main benefit of the module is that it allows Snakemake to submit all jobs in the pipeline immediately, rather than waiting for jobs to complete before releasing dependencies. This improves runtime dramatically when the SLURM queue is crowded. `run_pipeline.sh` handles submission of the pipeline to SLURM with sensible parameters.

## scripts/

Contains a bunch of single use analysis scripts used by the pipeline. Some of these scripts require my personal libraries ([fmbiopy](https://www.github.com/fennerm/fmbiopy) & [fen.R.util](https://www.github.com/fennerm/fen.R.util)) to be installed on the system.

## config/

This holds configuration files for external tools used in the pipeline.

## config.yml

This contains paths for relevant libraries and tools used by the pipeline.

## util/

This contains external binaries for tools which we didn't want to install globally on the user account. (Not included in the github repo due to licensing)

## input/ and output/
Input and output files of the pipeline.





