#!/usr/bin/Rscript

'Randomize the nucleotide order for each sequence in a fasta file

Usage: scramble_fasta.R [-o PATH] [--help] FASTA ...

-h --help    show this
-o --outdir=PATH    Output Directory (Default: Current directory)' -> doc


#' Scramble a fasta file
#'
#' For each sequence in the file, randomly scramble the sequence
#' @param fasta Input fasta file
#' @return A scrambled DNAstring object
#' @importFrom fen.R.util read_text ulapply
#' @importFrom Biostrings readDNAStringSet DNAStringSet
#' @export
scramble_fasta <- function(fasta) {
    fasta_string <- readDNAStringSet(fasta)
    DNAStringSet(lapply(fasta_string, sample))
}

# If run as script, scramble all command line arguments
if (!interactive()) {
    library(Biostrings)
    library(fen.R.util)
    library(tools)
    library(docopt)

    opts <- docopt(doc)

    # Set output directory
    if (is.null(opts$`--outdir`)) {
        outdir <- getwd()
    } else {
        outdir <- opts$`--outdir`
    }

    lapply(opts$FASTA, function(x) {

        outfile <- paste0(outdir, "/", file_path_sans_ext(x),
                          ".scrambled.fasta")
        scrambled <- scramble_fasta(x)
        writeXStringSet(scrambled, outfile)
    })
}
