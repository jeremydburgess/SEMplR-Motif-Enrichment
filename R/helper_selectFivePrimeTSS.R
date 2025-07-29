#' Select the 5′-most TSS per gene
#'
#' @param df A data.frame or tibble with columns `entrez`, `tss`, and `tx_strand`.
#' @return A tibble filtered to one row per `entrez`, picking the 5′-most TSS
#'   (respecting strand) for each gene.
#' @importFrom dplyr group_by arrange if_else desc slice_head ungroup
#' @keywords internal
helper_selectFivePrimeTSS <- function(df) {
  df %>%
    dplyr::group_by(entrez) %>%
    dplyr::arrange(
      entrez,
      dplyr::if_else(tx_strand == "+",
                     dplyr::desc(tss),  # for “+” strand pick highest (5′)
                     tss           # for “–” strand pick lowest (5′)
      )
    ) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()
}


