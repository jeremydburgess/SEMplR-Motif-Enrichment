#' Choose the Appropriate src_organism Table for a Given Keytype
#'
#' @description
#' \code{helper_chooseTable()} inspects all tables in a \code{src_organism} object
#' and returns the name of the table best suited for converting a set of IDs
#' of a given keytype (e.g., \dQuote{ensembltrans}, \dQuote{entrez}, etc.).
#' It applies priority rules to prefer transcript‐level tables when requested.
#'
#' @param so_obj     A \code{src_organism} object (from \pkg{Organism.dplyr}).
#' @param best       Character scalar, the keytype/column name you wish to map.
#' @param transcript Logical; if \code{TRUE}, prefer transcript‐level tables.
#'
#' @return A character string giving the name of the selected table within
#'   \code{so_obj} (e.g., \dQuote{id}, \dQuote{id_transcript}, etc.).
#'
#' @examples
#' \dontrun{
#' mapping <- buildMappingObject("Homo sapiens")
#' chooseTable(mapping$so_obj, best = "ensembltrans", transcript = TRUE)
#' }
#'
#' @keywords internal
#' @importFrom dplyr src_tbls tbl_vars
#' @export
helper_chooseTable <- function(so_obj, best, transcript) {
  tables <- src_tbls(so_obj)
  has_key <- vapply(tables, function(tbl_nm) {
    best %in% tbl_vars(tbl(so_obj, tbl_nm))
  }, logical(1))
  candidates <- tables[has_key]

  # priority rules:
  # 1) If our key is ensembltrans, force the transcript table
  if (best == "ensembltrans" && "id_transcript" %in% candidates) {
    return("id_transcript")
  }
  # 2) If we're in transcript mode, prefer transcript table
  if (transcript && "id_transcript" %in% candidates) {
    return("id_transcript")
  }
  # 3) Otherwise prefer the generic id (gene) table
  if ("id" %in% candidates) {
    return("id")
  }
  # 4) Fallback to the first table that contains your key
  if (length(candidates) > 0) {
    return(candidates[[1]])
  }
  stop("No table contains column '", best, "'")
}
