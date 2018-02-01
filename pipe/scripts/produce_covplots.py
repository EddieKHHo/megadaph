"""Snakescript for producing a series of covplots"""
from plumbum import local
from snakemake import shell

shell("blobtools covplot -p 20 -c {input.covplot} -i {input.db} -o {sample}")
outdir = local.path("output/produce_covplots1")
output_files = outdir.glob(snakemake.wildcards.sample + "*")
for f in output_files:
    extensions = ''.join(f.split(".")[2:])
    newname = outdir / (snakemake.wildcards.sample + extensions)
    f.move(newname)
