#' Randomly sample background regions
#'
#' @param pool A GRanges of candidate bg regions.
#' @param focal A GRanges of foreground regions.
#' @param n_ratio Numeric; number of bg = n_ratio * number of fg.
#' @param seed Integer; passed to set.seed() for reproducibility (default NULL)
#' @return A GRanges of randomly sampled background regions.
#' @export
helper_randomBackground <- function(pool,
                                    focal,
                                    n_ratio,
                                    bgReplace = FALSE,
                                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_fg   <- length(focal)
  n_pool <- length(pool)
  n_bg   <- n_ratio * n_fg

  if (n_bg > n_pool) {
    stop(
      "Requested ", n_bg, " background regions (",
      n_ratio, "×", n_fg, ") but only ",
      n_pool, " available in the pool."
    )
  }

  idx <- sample(
    seq_len(n_pool),
    size    = n_bg,
    replace = bgReplace
  )
  pool[idx]
}
