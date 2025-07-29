#' Map User IDs to Entrez and Determine Best Keytype
#'
#' Given a vector of user-supplied gene/transcript IDs, finds the AnnotationDbi
#' keytype (e.g. ens, refseq, symbol, etc.) that maps the highest fraction of
#' inputs, and returns both foreground and full background sets as Entrez IDs.
#' Optionally collapses transcript‐style inputs to genes when requested or when
#' reverse‐mapping inflation exceeds a threshold.
#'
#' @param foreground_ids Character vector of input IDs (e.g. ENSG, ENST, NM_*)
#' @param mapping        A list from \code{buildMappingObject()}, containing:
#'   \itemize{
#'     \item \code{so_obj}: the \code{src_organism} object
#'     \item \code{orgdb}:  the loaded OrgDb package object
#'     \item \code{organism}, \code{genomeBuild}, \code{txdb}: parameters used
#'   }
#' @param threshold      Fraction in [0,1]; minimum mapping rate to accept a
#'   keytype without falling back (default 0.9).
#' @param transcript     Logical; if \code{TRUE}, analyze as transcript‐level IDs
#'   (default \code{FALSE}).
#' @param stripVersions  Logical; strip trailing “.1”, “.2” from IDs (default \code{TRUE}).
#' @param inflateThresh  Fraction in [0,1]; if reverse‐mapping shows excessive inflation,
#'   automatically collapse transcripts to genes (default 1 ie. 100%).
#'
#' @return
#' A named list combining the original \code{mapping} components with:
#' \describe{
#'  \item{\code{fg_ids}}{data.frame(entrez, mappedID) for your foreground set}
#'  \item{\code{bg_ids}}{data.frame(entrez, mappedID) for the full background}
#'  \item{\code{userIDtype}}{the chosen keytype (e.g. "ensembl")}
#'  \item{\code{transcript}}{logical, whether transcript‐level mapping was used}
#' }
#'
#' @examples
#' \dontrun{
#'   mapping <- buildMappingObject("Homo sapiens")
#'   out     <- mapForegroundIDs(
#'     foreground_ids = c("ENSG00000139618","ENSG00000157764"),
#'     mapping        = mapping
#'   )
#'   str(out)
#' }
#'
#' @importFrom AnnotationDbi columns keytypes
#' @importFrom dplyr tbl filter select distinct collect group_by n_distinct
#' @importFrom rlang sym
#' @keywords internal
#' @export
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

  message("attempting to map foreground_ids to high likelihood columns...")
  fast_cols   <- intersect(cand_all, helper_guessCols(foreground_ids))

  stats_fast <- data.frame(
    id_type = fast_cols,
    matched = vapply(
      fast_cols,
      function(col) helper_countMapped(so_obj, foreground_ids, col),
      integer(1)
    ),
    stringsAsFactors = FALSE
  )

  # Calculate mapping percentages and convert NA -> 0
  stats_fast <- helper_cleanStats(stats_fast)

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
          function(col) helper_countMapped(so_obj, foreground_ids, col),
          integer(1)
        ),
        stringsAsFactors = FALSE
      )

    # Calculate mapping percentages and convert NA -> 0
    stats_slow <- helper_cleanStats(stats_slow)

    # Combine fast and slow stats results
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

  # Identify best table to use for converting ids
  tbl_nm <- helper_chooseTable(so_obj, best, transcript)

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
  # transcript = FALSE and warn they aey are losing transcript level coordinate specificity with this flag.
  #
  # Additionally reverse the mapping (entrez -> mapped id) to detect other transcript-
  # style ids by  1->many gene to mapped_id inflations and warn again.
  if (!transcript && best == "ensembltrans") {
    warning(
      "It looks like you provided Ensembl transcript IDs (", sQuote(best),
      ") but requested gene‐level analysis (transcript = FALSE).\n",
      "Any downstream coordinate lookup will use gene (Entrez) IDs, so you’ll lose\n",
      "the per‐transcript specificity of your input. If you really want transcript-\n",
      "level coordinates, set transcript = TRUE or supply Ensembl gene IDs instead."
    )

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
          "Your IDs appear transcript‐like: reverse‐mapping shows %.1f%% more\n",
          inflation * 100
        ),
        sprintf(
          "(%d unique transcript IDs for %d genes). Downstream, only gene‐level\n",
          nMapped, nGenes
        ),
        "coordinates will be used. If you need transcript‐level analyses, set\n",
        "transcript = TRUE and use Ensembl transcript IDs."
      )
      }
    }
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
