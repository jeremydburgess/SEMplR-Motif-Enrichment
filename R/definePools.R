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
definePools <- function(
    mapped_data,
    id_level    = c("auto","gene","transcript"),
    gene_filter = NULL,
    universe_filter = NULL,
    universe_filter_logic  = c("or","and")
){

  # 1 - Process the inputs and arguments

    # - normalize arguments
    id_level              <- match.arg(id_level)
    universe_filter_logic <- match.arg(universe_filter_logic)


    # figure out what id_level was detected during matchUserList()
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

    # - pick a sensible default gene_filter if none supplied
      if (is.null(gene_filter)) {
        gene_filter <- if (id_level == "gene") "gene" else "allTSS"
      }

    # - Extract mapped dataframe from input list object and ensure TSS has been computed
    univ <- mapped_data$input_df
    if (!"TSS" %in% names(univ)) {
      univ <- univ %>%
        mutate(TSS = if_else(strand == "+", start, end))
    }


  # 2) apply universe_filter if provided
  # Next filtering steps restricting the universe (gene_type) and transcripts/genes included (gene_filter) can be applied to the background only.
  # Can later harmonize the foreground to this population by using a semi-join
    if (!is.null(universe_filter)) {
      cols <- names(universe_filter)

      if (universe_filter_logic == "and") {
        for (col in cols) {
          vals <- universe_filter[[col]]
          if (length(vals) > 1) {
            # Sequentially require each value in turn
            for (v in vals) {
              univ <- univ %>% filter(.data[[col]] == v)
            }
          } else {
            # Single‐value case stays as %in%
            univ <- univ %>% filter(.data[[col]] %in% vals)
          }
        }
      }

      else {
        # OR logic: keep any row matching any of the filters
        if (length(cols) > 1) {
          warning(
            "Multiple universe_filter criteria supplied: ",
            "combining with OR by default."
          )
        }
        keep <- Reduce(
          `|`,
          lapply(cols, function(col) {
            vals <- universe_filter[[col]]
            univ[[col]] %in% vals
          })
        )
        univ <- univ[keep, ]
      }
    }

  # 3 - Set mode for gene_filter filtering logic (gene, flag, modalTSS, uniqueTSS, allTSS)
  # validate/parse gene_filter
  gf <- gene_filter

  # First check if flag column, ensure it describes a real flag, and identify which column it refers to
  if (grepl("^flag:", gf)) {
    # strip off the “flag:” prefix to get the suffix
    suffix   <- sub("^flag:", "", gf)
    # build the real column name
    flag_col <- paste0("flag_", suffix)

    # which flag_ columns actually exist?
    all_flags    <- names(univ)[startsWith(names(univ),"flag_")]
    avail_suffix <- sub("^flag_","", all_flags)

    # error if they asked for one you don't have
    if (!flag_col %in% names(univ)) {
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

  # 4 - Set filtering logic for each mode:
  # 'gene'
  if (mode == "gene") {
    if (id_level != "gene") {
      warning(
        "You requested gene_filter = 'gene' but id_level = '",
        id_level, "'.\n",
        "Filtering to feature_type = 'gene' anyway—this may return zero rows."
      )
    }
    out_univ <- univ %>%
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

    out_univ <- univ %>%
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
      univ_tr <- univ %>% filter(feature_type == "transcript")

      # 2) find the most frequent TSS per gene
      modal <- univ_tr %>%
        count(gene_id, TSS, name = "n") %>%
        group_by(gene_id) %>%
        slice_max(n, with_ties = FALSE) %>%
        ungroup() %>%
        dplyr::select(gene_id, TSS)

      # 3) join back to the full rows
      out_univ <- univ_tr %>%
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
    univ_tr <- univ %>% filter(feature_type == "transcript")

    # 2) collapse to one row per distinct (gene_id, TSS)
    out_univ <- univ_tr %>%
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
    out_univ <- univ %>% filter(feature_type == "transcript")
  } else {stop("Unknown gene_filter mode: ", mode)}

 # 5 - Harmonize mapped foreground with the newly defined background universe
  # - extract mapped dfs
  fg   <- mapped_data$mapped_df

  # semi join foreground with filtered background universe
  return(list(backgroundUniverse = out_univ,
            foregroundElements = fg %>% semi_join(out_univ, by="unique_id")))
}
