compute_match_stats <- function(user_list,
                                map_df,
                                idCols = NULL) {
  if (is.null(idCols)) {
    idCols <- c(
      "unique_id", "gene_name", "gene_id", "gene_id_base",
      "entrez_id", "HGNC_symbol", "HGNC_id", "transcript_name",
      "transcript_id", "transcript_id_base", "RefSeq_TxID", "RefSeq_PxID"
    )
  }
  # Only keep ID columns that actually exist
  cols <- intersect(idCols, names(map_df))
  total_ids <- length(user_list)

  stats <- lapply(cols, function(col) {
    # Unique non-NA values in this column
    vals <- unique(map_df[[col]])
    # Count how many user IDs matched
    n_matched <- sum(user_list %in% vals, na.rm = TRUE)
    # Compute percent
    pct <- n_matched / total_ids * 100

    data.frame(
      id_type      = col,
      matched      = n_matched,
      total        = total_ids,
      pct_matched  = pct,
      match_str    = paste0(n_matched, " of ", total_ids),
      stringsAsFactors = FALSE
    )
  })

  # Combine into one data.frame
  do.call(rbind, stats)
}
