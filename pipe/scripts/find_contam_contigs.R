#!/usr/bin/env Rscript

" Generate a list of contigs to filter from the initial assembly

This should include only the most obvious contamination candidates:
1. Contigs with high GC (> 0.6)
2. Contigs with bacterial and plant/algae taxonomic assignments.
3. Contigs with low coverage (< 10)

Usage:
  find_contam_contigs.R BLOBTABLE
" -> doc

import::from(dplyr, filter, select, "%>%")
source("lib/R/blobtools_util.R")

main <- function(blobtable) {
  exclusions <- filter(blobtable,
    superkingdom.t == "Bacteria" | grepl("ophyta", phylum.t) | GC > 0.6 |
      cov_sum < 10) %>%
    select(name)
  exclusions <- noquote(unlist(exclusions))
  exclusions_string <- paste(exclusions, collapse = "\n")
  write(exclusions_string, stdout())
}

if (!interactive()) {
  blobtable_file <- commandArgs(trailingOnly = TRUE)
  blobtable <- read_blobtable(blobtable_file)
  main(blobtable, output_filename)
}
