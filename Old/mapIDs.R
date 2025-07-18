mapIds_to_df <- function(
    user_list,
    species,
    stripVersion = TRUE,
    threshold    = 0.9,
    verbose      = TRUE
) {
  # (0) Strip Ensembl/RefSeq versions if needed
  if (stripVersion) {
    pat <- "^(ENS[GTPS]\\d+|N[MRP]_\\d+)\\.(\\d+)$"
    idx <- grepl(pat, user_list, ignore.case = TRUE)
    user_list[idx] <- sub(pat, "\\1", user_list[idx], ignore.case = TRUE)
  }

  # (1) Resolve species → OrgDb pkg & baseKey
  sp_in <- tolower(species)
  hits  <- vapply(.speciesMap, function(info) {
    any(sp_in == tolower(info$code),
        sp_in == tolower(info$common),
        sp_in == tolower(info$latin))
  }, logical(1))
  if (!any(hits))  stop("Unknown species: ", species)
  if (sum(hits)>1) stop("Ambiguous species: ", species)
  info    <- .speciesMap[[which(hits)]]
  orgdb   <- info$orgdb_pkg
  baseKey <- info$baseKey
  if (verbose) message("Using ", orgdb, " with keytype ", baseKey)

  # (2) Install/load OrgDb
  if (!requireNamespace(orgdb, quietly=TRUE)) {
    if (verbose) message("Installing ", orgdb, " …")
    BiocManager::install(orgdb, ask=FALSE)
  }
  requireNamespace(orgdb, quietly=TRUE)
  od <- get(orgdb, envir = asNamespace(orgdb))

  # (3) Candidate keytypes minus functional fields
  all_cols      <- AnnotationDbi::columns(od)
  valid_kts     <- AnnotationDbi::keytypes(od)
  non_id_fields <- c(
    "GO","GOALL","ONTOLOGY","ONTOLOGYALL","PATH","PMID",
    "PROSITE","EVIDENCE","EVIDENCEALL","PFAM","ENZYME",
    "MAP","OMIM","INTERPRO","COMMON","DESCRIPTION"
  )
  cand_all <- setdiff(intersect(all_cols, valid_kts), non_id_fields)

  # (4) Fast‐guess
  guessCols <- function(ul) {
    if (all(grepl("^[0-9]+$", ul)))        return("ENTREZID")
    if (all(grepl("^ENS(G|T|P)\\d+", ul))) return(c("ENSEMBL","ENSEMBLTRANS","ENSEMBLPROT"))
    if (all(grepl("^(N[MRP]_)", ul)))      return("REFSEQ")
    c("SYMBOL","ALIAS")
  }
  fast_cols   <- intersect(cand_all, guessCols(user_list))

  stats_fast <- data.frame(
    id_type = fast_cols,
    matched = vapply(fast_cols, function(col) {
      sum(!is.na(
        tryCatch(
          AnnotationDbi::select(od, keys=user_list, keytype=col, columns=col)[[col]],
          error = function(e) rep(NA, length(user_list))
        )
      ))
    }, integer(1)),
    stringsAsFactors = FALSE
  )
  stats_fast$pct_matched <- stats_fast$matched / length(user_list) * 100

  # (5) Fallback if needed
  if (max(stats_fast$pct_matched) >= threshold*100) {
    stats <- stats_fast
  } else {
    slow_cols  <- setdiff(cand_all, fast_cols)
    if (verbose) {
      message(
        "Fast‐guess only matched ",
        round(max(stats_fast$pct_matched),1),
        "% → trying ", length(slow_cols), " more keytypes…"
      )
    }
    stats_slow <- data.frame(
      id_type = slow_cols,
      matched = vapply(slow_cols, function(col) {
        sum(!is.na(
          tryCatch(
            AnnotationDbi::select(od, keys=user_list, keytype=col, columns=col)[[col]],
            error = function(e) rep(NA, length(user_list))
          )
        ))
      }, integer(1)),
      stringsAsFactors = FALSE
    )
    stats_slow$pct_matched <- stats_slow$matched / length(user_list) * 100
    stats <- rbind(stats_fast, stats_slow)
  }

  # (6) Pick best keytype
  best_pct <- max(stats$pct_matched)
  best     <- stats$id_type[which.max(stats$pct_matched)]
  if (best_pct < threshold*100) {
    stop("Unable to map ≥", threshold*100, "% of your IDs.")
  }
  if (verbose) {
    message(sprintf(
      "Chosen keytype '%s' (%.1f%% matched).",
      best, best_pct
    ))
  }

  # (7) Build mapping tibble with baseKey, mapped_id, GENETYPE, and mapped flag

  # (7a) universe of all baseKey values
  universe <- AnnotationDbi::keys(od, keytype = baseKey)
  df       <- tibble::tibble( !!baseKey := universe )

  # (7b) best‐ID mapping (all matches)
    mapped <- AnnotationDbi::select(
       od,
       keys    = user_list,
       keytype = best,
       columns = c(baseKey, best)
     ) %>%
      dplyr::rename(mapped_id = !!best)

  # (7c) optional GENETYPE
  if ("GENETYPE" %in% all_cols) {
    gt <- AnnotationDbi::select(
      od,
      keys      = universe,
      keytype   = baseKey,
      columns   = "GENETYPE",
      multiVals = "first"
    ) %>%
      dplyr::rename(!!baseKey := !!sym(baseKey))
    df <- dplyr::left_join(df, gt, by = baseKey)
  } else {
    df$GENETYPE <- NA_character_
  }

   # (7d) join in *all* mapped_id rows, producing one row per match
     df <- dplyr::left_join(df, mapped, by = baseKey) %>%
       dplyr::mutate(mapped = !is.na(mapped_id))

  df
}
