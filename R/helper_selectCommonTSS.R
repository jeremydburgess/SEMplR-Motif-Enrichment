#' Select the most common TSS per gene (without tie-breaking)
#'
#' @description
#' \code{helper_selectCommonTSS()} takes a data frame of transcript start
#' sites and returns, for each gene (`entrez`), **all** TSS values that
#' occur at the maximal frequency.  Tie-breaking (e.g. picking a single
#' 5′-most TSS) is intentionally **not** performed here and can be applied
#' separately downstream.
#'
#' @param df A data.frame or tibble with at least these columns:
#'   - \code{entrez}    : gene identifiers
#'   - \code{tss}       : numeric transcription start sites
#'   - \code{tx_strand} : strand (“+” or “–”)
#'
#' @return A tibble containing one or more rows per \code{entrez}, each
#'   TSS having the highest observed frequency for that gene.
#'
#' @details
#' This helper does:
#' 1. Count how many times each \code{tss} appears per \code{entrez}.
#' 2. Keep only those \code{tss} where \code{freq == max(freq)}.
#' 3. Re-join any additional columns (e.g.\ \code{tx_strand}, transcript IDs).
#'
#' It does *not* resolve ties when multiple TSS share the same maximal
#' frequency; that can be handled later (e.g.\ picking the 5′-most).
#'
#' @importFrom dplyr count group_by filter ungroup select left_join
#' @keywords internal
helper_selectCommonTSS <- function(df) {
  df <- df %>%
    dplyr::count(entrez, tss, name = "freq") %>%
    dplyr::group_by(entrez) %>%
    dplyr::filter(freq == max(freq)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-freq) %>%
    dplyr::left_join(df, by = c("entrez", "tss")) %>%
    # dplyr::group_by(entrez) %>%
    # # for + strand pick largest tss (5′); for – pick smallest
    # dplyr::arrange(
    #   entrez,
    #   dplyr::if_else(tx_strand == "+", dplyr::desc(tss), dplyr::asc(tss))) %>%
    # dplyr::slice_head(n = 1) %>%
    # dplyr::ungroup()

  return(df)
}
