import::from(fen.R.util, ulapply, read_table)
import::from(dplyr, filter)
import::from(stringr, regex, str_detect, str_which, str_replace)

#' Read a blobtable file into a data.frame
#' @export
read_blobtable <- function(blobtable_file) {
  blobtable <- read_table(
    blobtable_file, sep="\t", comment.char = "#", header = FALSE)
  colnames(blobtable) <- get_colnames(blobtable_file)
  blobtable
}


#' Strip a number from a character vector if it is preceded by a '.'
strip_dot_number <- function(x) {
  str_replace(x, regex("\\.[0-9]+"), "")
}


#' Get the column names from a blobtable file
get_colnames <- function(blobtable_file) {
  blobtable_text <- readLines(blobtable_file)
  colname_row <- blobtable_text[str_which(blobtable_text, "# name")]
  colname_row <- str_replace(colname_row, "# ", "")
  column_names <- unlist(strsplit(colname_row, "\t"))
  column_names <- strip_dot_number(column_names)
  column_names
}

