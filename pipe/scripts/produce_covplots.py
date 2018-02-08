"""Snakescript for producing a series of covplots"""
from plumbum import local
from snakemake import shell


command = ("blobtools covplot -p 20 -c {snakemake.input.cov} -i "
           "{snakemake.input.db} -r {snakemake.params.rank} "
           "-o {snakemake.params.output_prefix}")
if snakemake.params.exclude:
    command += " --exclude {snakemake.params.exclude}"

shell(command)

outdir = local.path(snakemake.params.output_prefix).dirname
output_files = outdir.glob(snakemake.wildcards.sample + "*")
for f in output_files:
    extensions = '.'.join(f.name.split(".")[-2:])
    newname = outdir / (snakemake.wildcards.sample + "." + extensions)
    f.move(newname)
