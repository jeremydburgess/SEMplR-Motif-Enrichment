mapForegroundIDs <- function(foreground_ids,
                             so_obj,
                             threshold = 0.9,
                             transcript = FALSE,
                             stripVersions = TRUE,
                             verbose = TRUE){

  # Strip Ensembl/RefSeq versions from foreground_ids if requested
  if (stripVersions) {
    pat <- "^(ENS[GTPS]\\d+|N[MRP]_\\d+)\\.(\\d+)$"
    idx <- grepl(pat, foreground_ids, ignore.case = TRUE)
    foreground_ids[idx] <- sub(pat, "\\1", foreground_ids[idx], ignore.case = TRUE)
  }

  # Identify potential mapping keys from so_obj
  all_cols      <- AnnotationDbi::columns(so_obj)
  valid_kts     <- AnnotationDbi::keytypes(so_obj)
  non_id_fields <- c(
    "go","goall","ontology","ontologyall","path","pmid",
    "prosite","evidence","evidenceall","pfam","enzyme",
    "map","omim","interpro","common","description"
  )
  cand_all <- setdiff(intersect(all_cols, valid_kts), non_id_fields)

  # Heuristically guess likely id columns to attempt mapping with based on foreground_ids format
  guessCols <- function(ul) {
    if (all(grepl("^[0-9]+$", ul)))        return("entrezid")
    if (all(grepl("^ENS(G|T|P)\\d+", ul))) return(c("ensembl","ensembltrans","ensemblprot"))
    if (all(grepl("^(N[MRP]_)", ul)))      return("refseq")
    c("symbol","alias")
  }


  fast_cols   <- intersect(cand_all, guessCols(foreground_ids))

  stats_fast <- data.frame(
    id_type = fast_cols,
    matched = vapply(fast_cols, function(col) {
      tryCatch(
        {
          # pull back the vector of mapped values
          vals <- AnnotationDbi::select(
            so_obj,
            keys    = foreground_ids,
            keytype = col,
            columns = col
          )
          # count *distinct* mapped IDs
          length(unique(vals))
        },
        error = function(e) NA_integer_
      )
    }, integer(1)),
    stringsAsFactors = FALSE
  )


  stats_fast$pct_matched <- stats_fast$matched / length(foreground_ids) * 100
  stats_fast$pct_matched[is.na(stats_fast$pct_matched)] <- 0
  stats_fast$matched[is.na(stats_fast$matched)] <- 0


  # If no successful mapping, try remaining columns
  if (max(stats_fast$pct_matched) >= threshold * 100) {
    stats <- stats_fast
  } else {
    slow_cols <- setdiff(cand_all, fast_cols)
    if (verbose) {
      message(
        "Fast‐guess only matched ",
        round(max(stats_fast$pct_matched), 1),
        "% → trying ", length(slow_cols), " more keytypes…"
      )
    }
    stats_slow <- data.frame(
      id_type = slow_cols,
      matched = vapply(
        slow_cols,
        function(col) {
          tryCatch(
            {
              vals <- AnnotationDbi::select(
                so_obj,
                keys    = foreground_ids,
                keytype = col,
                columns = col
              )
              length(unique(vals))
            },
            error = function(e) NA_integer_
          )
        },
        integer(1)  # one integer per keytype
      ),
      stringsAsFactors = FALSE
    )
    stats_slow$pct_matched <- stats_slow$matched / length(foreground_ids) * 100
    stats_slow$pct_matched[is.na(stats_slow$pct_matched)] <- 0
    stats_slow$matched[is.na(stats_slow$matched)] <- 0

    stats <- rbind(stats_fast, stats_slow)
    }

  # Pick best matching keytype (if any)
  best_pct <- max(stats$pct_matched)
  best     <- stats$id_type[which.max(stats$pct_matched)]
  if (best_pct < threshold*100) {
    stop("Unable to map ≥", threshold*100, "% of your IDs.")
  } else {
  if (verbose) {
    message(sprintf(
      "Chosen keytype '%s' (%.1f%% matched).",
      best, best_pct
    ))
  }
  }
          # Function to find the appropriate table to use for converting best mapped id
          chooseTable <- function(so_obj, best, transcript) {
            tables <- src_tbls(so_obj)
            has_key <- vapply(tables, function(tbl_nm) {
              best %in% tbl_vars(tbl(so_obj, tbl_nm))
            }, logical(1))
            candidates <- tables[has_key]

            # priority rules:
            if (transcript && "id_transcript" %in% candidates) {
              return("id_transcript")
            }
            if ("id" %in% candidates) {
              return("id")
            }
            if (length(candidates)) {
              return(candidates[1])
            }
            stop("No table contains column '", best, "'")
          }

  # Identify best table to use for converting ids
  tbl_nm <- chooseTable(so_obj, best, transcript)

    # Error is user is trying to use transcript mode on any id type except enstrans
    if(transcript && best != "ensembltrans"){
      stop(
        "Transcript-level analysis only supported with Ensembl transcript IDs.\n",
        "Your best keytype was ‘", best, "’. ",
        "Please convert your IDs to Ensembl transcript IDs (ensembltrans) and try again."
      )
    }

  # What if user specifies gene level (transcript = FALSE) analysis but provides transcript level ids (ensembltrans)?
  # Conceivable but rare. Warn and collapse to unique gene ids (entrez), retaining comma-separated
  # mapped transcripts associated with that gene id as mappedID
  if (!transcript && best == "ensembltrans") {

    warning(
      "You’ve mapped best to ", sQuote(best),
      " but requested gene-level analysis (transcript = FALSE).\n",
      "I will collapse transcripts → their parent genes for you, ",
      "but please double-check that this\nis what you intend."
    )

    # Pull and collapse transcript→gene for fg
    fg_ids <- tbl(so_obj, "id_transcript") %>%
      filter(ensembltrans %in% foreground_ids) %>%
      dplyr::select(entrez, mappedID = ensembltrans) %>%
      collect() %>%
      group_by(entrez) %>%
      summarise(
        mappedID = paste(unique(mappedID), collapse = ","),
        .groups = "drop"
      )

    # Pull and collapse transcript→gene for bg
    bg_ids <- tbl(so_obj, "id_transcript") %>%
      dplyr::select(entrez, mappedID = ensembltrans) %>%
      collect() %>%
      group_by(entrez) %>%
      summarise(
        mappedID = paste(unique(mappedID), collapse = ","),
        .groups = "drop"
      )

  } else {

    # Otherwise (ie gene level analysis from gene level id) pull entrez id and best
    # mapped id from appropriate table (1:1 mapping)

    # Foreground is only those records that map to foreground_ids
    fg_ids <- tbl(so_obj, tbl_nm) %>%
      filter(!!sym(best) %in% foreground_ids) %>%
      dplyr::select(entrez,mappedID = !!sym(best)) %>%
      distinct() %>% collect()

    #Background pool is all records from that table
    bg_ids <- tbl(so_obj, tbl_nm) %>%
      dplyr::select(entrez,mappedID = !!sym(best)) %>%
      distinct() %>% collect()

  }
  # Make output list
    mapped <- list(fg_ids = fg_ids,
                   bg_ids = bg_ids,
                   userIDtype = best,
                   transcript = transcript
                   )

  return(mapped)

  }
