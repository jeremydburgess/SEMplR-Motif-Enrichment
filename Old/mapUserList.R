source("R/computeMatchStats.R")

mapUserList <- function(
    user_list,
    map_df,
    idCols = NULL,
    threshold = 0.90) {

  # Calculate and tabulate mapping stats per id_type in map_df
  stats <- computeMatchStats(user_list, map_df, idCols)

  # Isolate stats for best mapping id_type
  best <- stats[which.max(stats$pct_matched), ]

  # If user_list doesn't map above defined threshold (default 90%) fail and request reformatted user_list
  if (best$pct_matched / 100 < threshold) {
    message("Unable to map IDs ≥ ", threshold*100, "% in any field; please supply IDs in a supported format.")
    return(NULL)
  }

  # Notify best mapped id and mapping statistics for it
  message(sprintf(
    "Identified input ids as %s (%d/%d; %.1f%%).",
    best$id_type, best$matched, best$total, best$pct_matched
  ))

  # Store best mapped id_type for output
  bestId <- best$id_type

  # Classify id_level (gene if in gene_key list. Otherwise, transcript)
  gene_keys <- c("gene_id","gene_name","entrez_id","HGNC_symbol","HGNC_id")
  # transcript_keys <- c("transcript_id","transcript_name","RefSeq_TxID") # Not actually used
  id_level <- if (bestId %in% gene_keys) "gene" else "transcript"


  # Ensure TSS has been created and create it if not
  df <- map_df

  df <- df %>%
    mutate(
      TSS = coalesce(
        TSS,                      # existing values (if any)
        if_else(strand == "+", start, end)  # fallback
      )
    )

  # Limit to minimum relevant columns
  # input_id, unique_id, transcript_id, gene_id, feature_type, gene_type, strand, flags, TSS, seqnames (chromosome)
  dfOut <- df %>%
    # 1) make mappedInput_id first
    mutate(
      mappedInput_id = !!sym(bestId)
    ) %>%
    # 2) then select everything you want
    dplyr::select(
      unique_id,
      mappedInput_id,
      gene_id,
      transcript_id,
      feature_type,
      gene_type,
      seqnames,
      strand,
      TSS,
      starts_with("flag_")
    )

  # Restrict map_df to only rows mapping to user list
  dfOut_mapped <- dfOut %>%
    filter(.data[[bestId]] %in% user_list)

  # return all three in a named list:
  list(
    match_stats   = stats,
    best_id_type  = bestId,
    id_level      = id_level,
    input_df     = dfOut,
    mapped_df     = dfOut_mapped
    )
  }

