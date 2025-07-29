#' Convert full organism name to BSgenome code
#'
#' Takes a two‐word Latin species name (e.g. “Homo sapiens”) and returns the
#' BSgenome convention of a single capitalized genus initial plus the species
#' name (e.g. “Hsapiens”).
#'
#' @param organism Character(1). A full species name in the form “Genus species”.
#'   Extra spaces between words are tolerated.
#' @return Character(1). A code composed of the uppercase first letter of the
#'   genus concatenated with the species epithet.
#' @keywords internal
helper_bsOrganismCode <- function(organism) {
  parts <- strsplit(organism, " +")[[1]]
  # Take first letter of Genus, then full Species
  code  <- paste0(substr(parts[1], 1, 1), parts[2])
  # ensure leading capital
  substr(code, 1, 1) <- toupper(substr(code, 1, 1))
  code
}

