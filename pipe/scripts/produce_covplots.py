"""Snakescript for producing a series of covplots"""
from plumbum import local
from snakemake import shell

shell("blobtools covplot -p 20 -c {snakemake.input.covplot} -i "
      "{snakemake.input.db} -o {snakemake.params.output_prefix}")
outdir = local.path("output/produce_covplots1")
output_files = outdir.glob(snakemake.wildcards.sample + "*")
for f in output_files:
    extensions = '.'.join(f.name.split(".")[2:])
    newname = outdir / (snakemake.wildcards.sample + "." + extensions)
    f.move(newname)
