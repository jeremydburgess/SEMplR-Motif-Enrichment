mapForegroundIDs <- function(foreground_ids,
                             mapping,
                             threshold     = 0.9,
                             transcript    = FALSE,
                             stripVersions = TRUE,
                             inflateThresh = 1
                             ) {

  # Pull out so_obj and and orgdb from mapping input list
  so_obj <- mapping$so_obj
  orgdb  <- mapping$orgdb

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

  message("attempting to map foreground_ids to high likelihood columns...")
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

      message(
        "Fast‐guess only matched ",
        round(max(stats_fast$pct_matched), 1),
        "% → trying ", length(slow_cols), " more keytypes…"
      )

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

    message(sprintf(
      "Chosen keytype '%s' (%.1f%% matched).",
      best, best_pct
    ))

  }
          # Function to find the appropriate table to use for converting best mapped id
          chooseTable <- function(so_obj, best, transcript) {
            tables <- src_tbls(so_obj)
            has_key <- vapply(tables, function(tbl_nm) {
              best %in% tbl_vars(tbl(so_obj, tbl_nm))
            }, logical(1))
            candidates <- tables[has_key]

            # priority rules:
            # 1) If our key is ensembltrans, force the transcript table
            if (best == "ensembltrans" && "id_transcript" %in% candidates) {
              return("id_transcript")
            }
            # 2) If we're in transcript mode, prefer transcript table
            if (transcript && "id_transcript" %in% candidates) {
              return("id_transcript")
            }
            # 3) Otherwise prefer the generic id (gene) table
            if ("id" %in% candidates) {
              return("id")
            }
            # 4) Fallback to the first table that contains your key
            if (length(candidates) > 0) {
              return(candidates[[1]])
            }
            stop("No table contains column '", best, "'")
          }

  # Identify best table to use for converting ids
  tbl_nm <- chooseTable(so_obj, best, transcript)

    # Error if user is trying to use transcript mode on any id type except ensembltrans
    if(transcript && best != "ensembltrans"){
      stop(
        "Transcript-level analysis only supported with Ensembl transcript IDs.\n",
        "Your best keytype was ‘", best, "’. ",
        "Please convert your IDs to Ensembl transcript IDs (ensembltrans) and try again."
      )
    }

  # Otherwise pull entrez id and best mapped id from appropriate table for both
  # mapped ids (fg), and entire table population (bg)

  # Foreground is only those records that map to foreground_ids
  fg_ids <- tbl(so_obj, tbl_nm) %>%
    dplyr::filter(!!sym(best) %in% foreground_ids) %>%
    dplyr::select(entrez,mappedID = !!sym(best)) %>%
    dplyr::distinct() %>%
    dplyr::collect()

  #Background pool is all records from that table
  bg_ids <- tbl(so_obj, tbl_nm) %>%
    dplyr::select(entrez,mappedID = !!sym(best)) %>%
    dplyr::distinct() %>%
    dplyr::collect()

  # Ensure background pool is the same universe as foreground by dropping any rows
  # with no available value for best mappedID type
  bg_ids <- bg_ids %>%
    dplyr::filter(!is.na(mappedID))

  # What if user specified gene level analysis (transcript = FALSE) but
  # foreground_ids provided are for transcripts.
  # Conceivable but rare.
  #
  # Specifically detect most common version of that case (best = ensembltrans &
  # transcript = FALSE and collapse by gene_id (entrez) (with warning).
  #
  # Additionally reverse the mapping (entrez -> mapped id) to detect other transcript-
  # style ids by  1->many gene to mapped_id inflations and collapse these too (with warning).

  doCollapse = FALSE
  # test for ensembltrans + transcript = FALSE
  if (!transcript && best == "ensembltrans") {
    warning(
      "You’ve mapped best to ", sQuote(best),
      " but requested gene-level analysis (transcript = FALSE).\n",
      "I will collapse transcripts → their parent genes for you, ",
      "but please double-check that this\nis what you intended."
    )
    doCollapse = TRUE
  } else if (!transcript){

    # test for other instances of 1:many gene:foreground_id inflations
    # reverse the mapping from the fg data
    rev_df <- tbl(so_obj, tbl_nm) %>%
      dplyr::filter(entrez %in% fg_ids$entrez) %>%     # raw user IDs not yet collapsed
      dplyr::select(entrez, mappedID = !!sym(best)) %>%
      dplyr::collect()

    # count genes, and mapped_ids and look for excessive inflation
    nGenes   <- n_distinct(rev_df$entrez)
    nMapped  <- n_distinct(rev_df$mappedID)
    if (nGenes > 0) {
      inflation <- nMapped / nGenes - 1
    if (inflation > inflateThresh) {
      warning(
        sprintf(
          paste0(
            "You have requested gene-level analysis but it looks like the IDs provided\n",
            "represent transcripts. Reverse mapping from gene ids back to your provided\n",
            "ids inflates IDs by %.1f%% (%d mappedIDs for %d genes).\n\n",
            "I will collapse transcripts → their parent genes for you, but please\n",
            "double-check that this is what you intended."
          ),
          inflation * 100,
          nMapped,
          nGenes
        ))
      doCollapse <- TRUE
      }
    }
  }

  # 4) If triggered, collapse both fg and bg by entrezid, retaining all mappedIDs
  # as a coma-separaated string in mappedID
  if (doCollapse) {
    fg_ids <- fg_ids %>%
      group_by(entrez) %>%
      summarise(
        mappedID = paste(unique(mappedID), collapse = ","),
        .groups  = "drop"
      )
    bg_ids <- bg_ids %>%
      group_by(entrez) %>%
      summarise(
        mappedID = paste(unique(mappedID), collapse = ","),
        .groups  = "drop"
      )
  }

  # Make output list
  mapped <- c(
    mapping,         # everything that came in (so_obj, orgdb, organism, genomeBuild, txdb, etc.)
    list(
      fg_ids     = fg_ids,
      bg_ids     = bg_ids,
      userIDtype = best,
      transcript = transcript
    )
  )

  return(mapped)

  }
