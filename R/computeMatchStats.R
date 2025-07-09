compute_match_stats <- function(user_list,
                                map_df,
                                idCols = NULL) {
  total_ids <- length(user_list)

  # 1) Pick candidate columns
  if (is.null(idCols)) {
    idCols <- names(map_df)[
      vapply(map_df, function(col) {
        is.character(col) || is.factor(col) || is.integer(col) || is.numeric(col)
      }, logical(1))
    ]
  }
  cols <- intersect(idCols, names(map_df))

  # 2) Loop over each potential ID column
  stats <- lapply(cols, function(col) {
    vals <- unique(map_df[[col]])

    # if vals are numeric, coerce user_list to numeric for matching
    if (is.numeric(vals) || is.integer(vals)) {
      user_ids_num <- suppressWarnings(as.numeric(user_list))
      n_matched   <- sum(user_ids_num %in% vals, na.rm = TRUE)
    } else {
      n_matched   <- sum(user_list %in% vals, na.rm = TRUE)
    }

    pct <- n_matched / total_ids * 100

    data.frame(
      id_type     = col,
      matched     = n_matched,
      total       = total_ids,
      pct_matched = pct,
      match_str   = sprintf("%d of %d", n_matched, total_ids),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, stats)
}
