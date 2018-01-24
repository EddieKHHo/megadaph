
######## Snakemake header ########
import sys; sys.path.insert(0, "/mnt/lfs2/schaack/.linuxbrew/opt/python3/lib/python3.6/site-packages/snakemake-4.3.1-py3.6.egg"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05X"\x00\x00\x00output/initial_assembly/GASC.fastaq\x06a}q\x07X\x06\x00\x00\x00_namesq\x08}q\tsbX\x06\x00\x00\x00outputq\ncsnakemake.io\nOutputFiles\nq\x0b)\x81q\x0cX\x1e\x00\x00\x00output/haplomerger2/GASC.fastaq\ra}q\x0eh\x08}q\x0fsbX\x06\x00\x00\x00paramsq\x10csnakemake.io\nParams\nq\x11)\x81q\x12X\x02\x00\x00\x0010q\x13a}q\x14(h\x08}q\x15X\x07\x00\x00\x00threadsq\x16K\x00N\x86q\x17sh\x16h\x13ubX\t\x00\x00\x00wildcardsq\x18csnakemake.io\nWildcards\nq\x19)\x81q\x1aX\x04\x00\x00\x00GASCq\x1ba}q\x1c(h\x08}q\x1dX\x06\x00\x00\x00sampleq\x1eK\x00N\x86q\x1fsX\x06\x00\x00\x00sampleq h\x1bubh\x16K\x01X\t\x00\x00\x00resourcesq!csnakemake.io\nResources\nq")\x81q#(K\x01K\x01e}q$(h\x08}q%(X\x06\x00\x00\x00_coresq&K\x00N\x86q\'X\x06\x00\x00\x00_nodesq(K\x01N\x86q)uh&K\x01h(K\x01ubX\x03\x00\x00\x00logq*csnakemake.io\nLog\nq+)\x81q,X \x00\x00\x00output/haplomerger2/log/GASC.logq-a}q.h\x08}q/sbX\x06\x00\x00\x00configq0}q1(X\x08\x00\x00\x00isolatesq2X\x03\x00\x00\x00allq3X\x08\x00\x00\x00readsdirq4X\x0b\x00\x00\x00input/readsq5X\x08\x00\x00\x00adaptersq6X#\x00\x00\x00input/adapters/adapters_illumina.faq7X\n\x00\x00\x00threadsmaxq8K\nX\x11\x00\x00\x00spades_executableq9X\x0f\x00\x00\x00~/bin/spades.pyq:uX\x04\x00\x00\x00ruleq;X\x0c\x00\x00\x00haplomerger2q<ub.'); from snakemake.logging import logger; logger.printshellcmds = True
######## Original script #########
#!/usr/bin/env python
"""Snakemake script for running Haplomerger2"""
import sys

from plumbum import (
    FG,
    local,
)
from plumbum.cmd import (
    gzip,
    windowmasker,
)
from snakemake import shell

assembly = local.path(snakemake.input[0])
try:
    threads = snakemake.params.threads
except:
    threads = "1"

if not snakemake.wildcards.sample:
    sys.exit(1)

import pdb; pdb.set_trace()
# Copy the template project directory
template = local.path("lib", "haplomerger2", "template")
projdir = local.path("lib", "haplomerger2", snakemake.wildcards.sample)
sys.path.append(str(projdir))
if projdir.exists():
    projdir.delete()
template.copy(projdir)

# Create a masking library
masklib = projdir / "masklib.ustat"
windowmasker["-checkdup", "true", "-mk_counts", "-in", assembly, "-out",
             masklib, "-mem", "10000"] & FG

# Soft mask assembly
masked_assembly = projdir / (snakemake.wildcards.sample + ".fa")
windowmasker["-ustat", masklib, "-in", assembly, "-out", masked_assembly,
             "-outfmt", "fasta", "-dust", "true"] & FG
gzip[masked_assembly] & FG
masked_assembly = local.path(masked_assembly + ".gz")

with local.cwd(projdir):
    shell(' '.join(["./run_all.batch", snakemake.wildcards.sample, threads,
                    "2>", snakemake.log]))
