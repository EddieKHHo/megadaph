#!/usr/bin/env Rscript
" Generate a list of contigs to filter from the second assembly iteration

Usage:
  find_contam_contigs.R BLOBTABLE
" -> doc

# This repeats the the targets from the first filter, namely:
# 1. Contigs with high GC (> 0.6)
# 2. Contigs with bacterial and plant/algae taxonomic assignments.
# 3. Contigs with low coverage (< 10)
#
# Some of the above criteria were refined slightly from the first filtration
# step.
#
# In addition the following are also targeted:
# 1. Contigs with coverage < 5 in any single sample
# 2. Contigs from fungal all fungal phyla except microsporidia


library(dplyr)
source("lib/R/blobtools_util.R")

main <- function(blobtable) {
  exclusions1 <- filter(blobtable,
    superkingdom.t %in% c("Bacteria", "Viruses", "Archaea") |
      grepl("ophyta", phylum.t) |
      phylum.t == "Eustigmatophyceae" |
      grepl("mycota", phylum.t) |
      GC > 0.6 |
      cov_sum < 10 |
      order.t == "Primates" |
      order.t == "Rodentia" |
      order.t == "Carnivora") %>%
    select(name)
  exclusions2 <- filter_at(blobtable,
    vars(starts_with("cov")), any_vars(. < 5)) %>%
    select(name)
  exclusions <- c(noquote(unlist(exclusions1)), noquote(unlist(exclusions2)))
  exclusions <- unique(exclusions)
  exclusions_string <- paste(exclusions, collapse = "\n")
  write(exclusions_string, stdout())
}

if (!interactive()) {
  blobtable_file <- commandArgs(trailingOnly = TRUE)
  blobtable <- read_blobtable(blobtable_file)
  main(blobtable)
}
