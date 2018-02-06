#!/usr/bin/env Rscript

" Filter a produced from the first round of assemblies

This should remove the most obvious contamination candidates:
1. Contigs with high GC (> 0.6)
2. Contigs with bacterial and plant/algae taxonomic assignments.

Usage:
  filter_blobtable1.R --output=FILTERED_TABLE BLOBTABLE

Options:
  -o, --output=FILTERED_TABLE   Filename of the output table
" -> doc

import::from(docopt, docopt)
import::from(dplyr, filter, "%>%")
import::from(fen.R.util, write_table)
source("lib/R/blobtools_util.R")

main <- function(blobtable, output_filename) {
  filtered_table <- filter(blobtable, superkingdom.t != "Bacteria") %>%
    filter(GC < 0.6) %>%
    filter(cov_sum > 15)
  write_table(filtered_table, output_filename, sep = "\t")
}

if (!interactive()) {
  opts <- docopt(doc)
  blobtable_file <- unlist(opts["BLOBTABLE"])
  output_filename <- unlist(opts["output"])
  blobtable <- read_blobtable(blobtable_file)
  main(blobtable, output_filename)
}
