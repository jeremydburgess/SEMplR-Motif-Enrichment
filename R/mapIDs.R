#' @title Build an ID‐mapping table from aspecies-specific OrgDb package
#'
#' @description In response to the user-provided species, install and load the
#' relevant OrgDb package, mapping a base keytype (In most cases ENTREZID) to
#' all other #’ supported identifier columns.
#'
#' @param user_list A vector of identifiers to form the foreground for enrichment
#' analysis
#' @param species string defining species (accept abbreviation or full name)
#' @param threshold Minimal succesfful mapping proportion to avoid exiting the function
#'
#' @return A tibble with user_ids_, keytype id, mapped (boolean), GENETYPE
#' @export
mapIds <- function(
    user_list,
    species,
    threshold = 0.9) {



}
