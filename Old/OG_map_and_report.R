# OG Map and report for safe keeping


# map_and_report: Choose the best‐matching ID type for a user list and return minimal coordinates
#
# Args:
#   user_list:  character vector of gene or transcript IDs supplied by the user
#   map_df:     tibble produced by annotateGTF(), containing parsed GTF columns, flags, and external IDs
#   idCols:     optional character vector of column names in map_df to consider as possible ID matches;
#               if NULL, defaults to the standard set (gene_name, gene_id, transcript_id, etc.)
#   threshold:  numeric between 0 and 1 giving the minimum fraction of user_list that must match
#               (default 0.90) in the best‐hit column or the function returns NULL
#
# Returns:
#   If a match ≥ threshold is found: a tibble with columns
#     • unique_id, seqnames, start, end, strand
#     for only those rows whose ID (in the chosen column) is in user_list.
#   Otherwise: NULL (and a message advising the user to supply a different ID format).
map_and_report <- function(user_list, map_df, idCols = NULL, threshold = 0.90) {
  # 1) Compute match stats
  stats <- compute_match_stats(user_list, map_df, idCols)
  # stats <- compute_match_stats(user_list, gtf_map, idCols)

  # 2) Pick the best‐matching ID type
  best <- stats[which.max(stats$pct_matched), ]

  # 3) Check threshold
  if (best$pct_matched / 100 >= threshold) {
    msg <- sprintf(
      "Identified input ids as %s (%d of %d; %.1f%%).",
      best$id_type, best$matched, best$total, best$pct_matched
    )
    message(msg)

    # 4) Classify the ID type
    gene_keys       <- c("gene_id", "gene_name", "entrez_id", "HGNC_symbol", "HGNC_id")
    transcript_keys <- c("transcript_id", "transcript_name", "RefSeq_TxID")

    id_level <- if (best$id_type %in% gene_keys) {
      "gene"
    } else if (best$id_type %in% transcript_keys) {
      "transcript"
    } else {
      # fallback: assume transcript-level for anything else
      "transcript"
    }

    # 5) Subset to only that feature_type
    subset_df <- map_df %>%
      filter(feature_type == id_level,
             .data[[best$id_type]] %in% user_list)

    # 6) Build your minimal table
    minimal <- subset_df %>%
      select(
        unique_id, seqnames,
        start, end, strand
      ) %>%
      unique()

    return(minimal)
  } else {
    message(
      sprintf(
        "Unable to map IDs ≥ %.0f%% in any supported field; please supply IDs in one of: %s.",
        threshold*100,
        paste(stats$id_type, collapse = ", ")
      )
    )
    return(NULL)
  }
}
