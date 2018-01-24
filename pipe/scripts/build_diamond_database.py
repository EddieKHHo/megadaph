#!/usr/bin/env python
"""Script for downloading and initializing a Diamond database for Uniprot
reference proteomes

Usage:
  download_diamond_database.py OUTPUT_PREFIX
"""
import os

from plumbum import local
from snakemake import shell


def main(outdir, threads):
    with local.cwd(outdir):
        # Download database
        # shell("curl ftp://ftp.uniprot.org/pub/databases/uniprot/"
        #       "current_release/knowledgebase/reference_proteomes/"
        #       "Reference_Proteomes_2017_12.tar.gz > "
        #       "Reference_Proteomes_2017_12")

        # shell("tar -xvf Reference_Proteomes_2017_12")

        # # Unpack protein FASTAs for each kingdom
        # shell("ls */*.fasta.gz | grep -v 'DNA' | grep -v 'additional' | "
        #       "parallel -j" + threads + " gunzip {{}}")

        shell("cat */*.fasta > uniprot_ref_proteomes.fasta")

        # # Simplify sequence IDs
        shell("cat uniprot_ref_proteomes.fasta | sed -r "
              "'s/(^>sp\|)|(^>tr\|)/>/g' | cut -f1 -d\"|\" > temp")
        shell("mv temp uniprot_ref_proteomes.fasta")

        # Make Diamond database
        shell("diamond makedb --in uniprot_ref_proteomes.fasta -d "
              "uniprot_ref_proteomes")

        # Subset mapping file to only contain NCBI TaxID entries

        shell("ls */*.idmapping.gz | parallel -j" + threads + " gunzip {{}}")
        shell("cat */*.idmapping | grep 'NCBI_TaxID' > "
              "uniprot_ref_proteomes.taxids")


if __name__ == "__main__":
    outdir = os.path.dirname(snakemake.output[0])
    main(outdir, str(snakemake.threads))
