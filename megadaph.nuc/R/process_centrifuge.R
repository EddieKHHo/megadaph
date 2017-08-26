#' Convert taxonomic ID to taxonomic name
#' @param tax_id Centrifuge taxonomic ID
#' @param cent_report Centrifuge report output
#' @return Taxonomic name
#' @export
tax_id_to_name <- function(tax_id, cent_report) {
    if (tax_id == 0) {
        "unknown"
    } else {
        cent_report[which(cent_report$taxID == tax_id), "name"]
    }
}

#' Add taxonomic names to centrifuge table
#' @param cent_out Centrifuge table
#' @param cent_report Centrifuge report table
#' @return Centrifuge table with names column
#' @importFrom fen.R.util read.tsv ulapply
#' @export
add_names_col <- function(cent_out, cent_report) {
    tax_id <- cent_out$taxID
    names <- ulapply(tax_id, tax_id_to_name, cent_report)
    cbind(cent_out, names)
}


add_match_fraction_col <- function(cent_out) {
    match_fraction <- cent_out$hitLength / cent_out$queryLength
    cbind(cent_out, match_fraction)
}

cov_from_read_id <- function(read_id) {
    underscore_split <- unlist(strsplit(read_id, "_"))
    coverage <- underscore_split[length(underscore_split)]
    as.numeric(coverage)
}

add_coverage_col <- function(cent_out) {
    coverage <- ulapply(cent_out$readID, cov_from_read_id)
    cbind(cent_out, coverage)
}

add_out_cols <- function(cent_out) {


}
