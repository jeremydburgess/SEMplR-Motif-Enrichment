#' Count distinct, non-NA mappings for a given keytype
#'
#' @description
#' Given a \code{src_organism} object and a vector of keys, attempts to
#' map each key through \code{AnnotationDbi::select()} for the specified
#' \code{keytype}, then counts how many unique, non-NA values were returned.
#'
#' @param so_obj    A \code{src_organism} data source (from Organism.dplyr).
#' @param keys      Character vector of input IDs to map (e.g. gene or transcript IDs).
#' @param keytype   Character scalar specifying the AnnotationDbi keytype/column.
#'
#' @return
#' An integer giving the number of distinct, non-NA mappings.
#'
#' @keywords internal
helper_countMapped <- function(so_obj, keys, keytype) {
  cnt <- tryCatch({
    vals <- AnnotationDbi::select(
      so_obj,
      keys    = keys,
      keytype = keytype,
      columns = keytype
    )
    length(unique(vals))
  }, error = function(e) {
    NA_integer_
  })
  return(cnt)
}
