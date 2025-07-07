#’ Assemble a foreground set of TSSs for motif analysis
#’
#’ Choose how to collapse or select TSS coordinates per gene or transcript.
#’
#’ @param mapped_data A list returned by `mapUserList()`, with elements
#’   - `match_stats`: data.frame of mapping percentages
#’   - `best_id_type`: the column name that best matched
#’   - `id_level`: “gene” or “transcript”
#’   - `mapped_df`: the data.frame of matched rows
#’ @param id_level   One of `"auto"`, `"gene"`, or `"transcript"`.
#’   `"auto"` (the default) will use `mapped_data$id_level`.
#’ @param gene_filter Either
#’   - `"gene"`, `"modalTSS"`, `"uniqueTSS"`, `"allTSS"`, or
#’   - `"flag:<FLAG_NAME>"` to pick only transcripts with a given `flag_<FLAG_NAME> == 1`.
#’   Defaults to `"gene"` for gene‐level IDs or `"allTSS"` for transcript‐level IDs.
#’
#’ @return A tibble with columns
#’   - `unique_id`, `seqnames`, `TSS`, `strand`, plus any flags if used.
#’   Exactly one or more rows per gene/transcript depending on `gene_filter`.
#’ @export
assembleTSSforeground <- function(
    mapped_data,
    id_level    = c("auto","gene","transcript"),
    gene_filter = NULL){

  # 1 - Process id_level
  # normalize the choice
  id_level <- match.arg(id_level)

  # figure out what was detected
  detected <- mapped_data$id_level

  # if user overrode
  if (id_level != "auto") {
    if (id_level != detected) {
      warning(
        "You requested id_level = '", id_level,
        "' but data looks like '", detected,
        "'. Proceeding as '", id_level, "'."
      )
    }
  } else {
    # auto → respect the detected level
    id_level <- detected
    message("Processing as '", id_level,
    "' level inputs.")
  }

  # now id_level is exactly "gene" or "transcript"

  # 2 - Pick a sensible default gene_filter if none supplied
  if (is.null(gene_filter)) {
    gene_filter <- if (id_level == "gene") "gene" else "allTSS"
  }


  # 3 - Extract mapped dataframe from input list object and ensure TSS has been computed
  mapped_df <- mapped_data$mapped_df
  if (!"TSS" %in% names(mapped_df)) {
    mapped_df <- mapped_df %>%
      mutate(TSS = if_else(strand == "+", start, end))
  }

  # 4 - Set mode for filtering logic (gene, flag, modalTSS, uniqueTSS, allTSS)
  # validate/parse gene_filter
  gf <- gene_filter

  # First check if flag column, ensure it describes a real flag, and identify which column it refers to
  if (grepl("^flag:", gf)) {
    # strip off the “flag:” prefix to get the suffix
    suffix   <- sub("^flag:", "", gf)
    # build the real column name
    flag_col <- paste0("flag_", suffix)

    # which flag_ columns actually exist?
    all_flags    <- names(mapped_df)[startsWith(names(mapped_df),"flag_")]
    avail_suffix <- sub("^flag_","", all_flags)

    # error if they asked for one you don't have
    if (!flag_col %in% names(mapped_df)) {
      stop(
        "assembleTSSforeground(): unknown flag '", suffix, "'.\n",
        "Available flags: ", paste(avail_suffix, collapse = ", ")
      )
    }

    mode <- "flag"
    # Otherwise set to relevant alternate mode
  } else {
    # non-flag choices
    allowed <- c("gene","modalTSS","uniqueTSS","allTSS")
    mode    <- match.arg(gf, choices = allowed)
    flag_col <- NULL
  }

  # 5 - Set filtering logic for each mode:
  # 'gene'
  if (mode == "gene") {
    if (id_level != "gene") {
      warning(
        "You requested gene_filter = 'gene' but id_level = '",
        id_level, "'.\n",
        "Filtering to feature_type = 'gene' anyway—this may return zero rows."
      )
    }
    out_df <- mapped_df %>%
      filter(feature_type == "gene")
  }

  # 'flag'
  else if (mode == "flag") {
    # Warn if the user supplied transcript IDs,
    # since this strategy is really gene-centric (picking transcripts for each gene).
    if (id_level == "transcript") {
      warning(
        "You requested a gene_filter = 'flag:<FLAG>' strategy, ",
        "but id_level = 'transcript'.\n",
        "Filtering to flagged transcripts only (feature_type == 'transcript')."
      )
    }

    out_df <- mapped_df %>%
      filter(
        feature_type == "transcript",
        !!sym(flag_col) == 1L
      )
  }

  # 'modalTSS'
  else if (mode == "modalTSS") {
    if (id_level != "gene") {
        warning(
          "You requested gene_filter = 'modalTSS' but id_level = '",
          id_level, "'.\n",
          "Computing modal TSS per gene_id across all matched transcripts."
        )
      }

      # 1) only transcripts
      df_tr <- mapped_df %>% filter(feature_type == "transcript")

      # 2) find the most frequent TSS per gene
      modal <- df_tr %>%
        count(gene_id, TSS, name = "n") %>%
        group_by(gene_id) %>%
        slice_max(n, with_ties = FALSE) %>%
        ungroup() %>%
        select(gene_id, TSS)

      # 3) join back to the full rows
      out_df <- df_tr %>%
        inner_join(modal, by = c("gene_id","TSS")) %>%
        # Ensure we now only have one TSS per gene_id
        distinct(gene_id, .keep_all = TRUE)
  }

  # 'uniqueTSS'
  else if (mode == "uniqueTSS") {
    if (id_level != "gene") {
      warning(
        "You requested gene_filter = 'uniqueTSS' but id_level = '",
        id_level, "'.\n",
        "Finding unique TSSs per gene_id across all matched transcripts."
      )
    }

    # 1) only transcripts
    df_tr <- mapped_df %>% filter(feature_type == "transcript")

    # 2) collapse to one row per distinct (gene_id, TSS)
    out_df <- df_tr %>%
      distinct(gene_id, TSS, .keep_all = TRUE)
  }

  # 'allTSS'
  else if (mode == "allTSS") {
    if (id_level != "gene") {
      warning(
        "You requested gene_filter = 'allTSS' but id_level = '",
        id_level, "'.\n",
        "Finding all TSSs per gene_id across all matched transcripts."
      )
    }

    # 1) only transcripts - no further filtering needed - returning all TSS for mapped transcripts, including duplicates
    out_df <- mapped_df %>% filter(feature_type == "transcript")
  } else {stop("Unknown gene_filter mode: ", mode)}

return(out_df)
}
