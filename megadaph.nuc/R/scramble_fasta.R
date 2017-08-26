#!/usr/bin/Rscript

#' Scramble a fasta file
#'
#' For each sequence in the file, randomly scramble the sequence and write to
#' file.
#' @param fasta Input fasta file
#' @param outfile Output file path
#' @importFrom fen.R.util read_text ulapply
#' @importFrom Biostrings readDNAStringSet DNAStringSet writeXStringSet
#' @export
scramble_fasta <- function(fasta, outfile) {
    fasta_string <- readDNAStringSet(fasta)
    shuffled <- DNAStringSet(lapply(fasta_string, sample))
    writeXStringSet(shuffled, outfile)
}

# If run as script, scramble all command line arguments
if (!interactive()) {
    library(Biostrings)
    library(fen.R.util)
    library(tools)
    args = commandArgs(trailingOnly = TRUE)
    outdir <- "scrambled"
    dir.create(outdir, showWarnings = FALSE)
    lapply(args, function(x) {
        outfile <- paste0(file_path_sans_ext(x), ".scrambled.fasta")
        scramble_fasta(x, outfile)
    })
}
