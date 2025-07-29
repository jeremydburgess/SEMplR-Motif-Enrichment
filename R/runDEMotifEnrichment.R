runDEMotifEnrichment <- function(foreground_ids,
                               organism      = "Homo sapiens",
                               genomeBuild   = "auto",
                               txdb          = "auto",
                               transcript    = FALSE,
                               geneType      = NULL,
                               TSS.method = c(
                                 "UCSCgene",
                                 "Ensembl_canonical",
                                 "commonTSS",
                                 "uniqueTSS",
                                 "fivePrimeTSS",
                                 "allTSS"
                               ),
                               threshold     = 0.9,
                               stripVersions = TRUE,
                               inflateThresh = 1,
                               standardChrom = TRUE,
                               overlapReduction = TRUE,
                               overlapMinGap = 0,
                               onePromoterPerGene = TRUE



) {

  TSS.method <- match.arg(TSS.method)

  # 1) Build the mapping object (cheap metadata + package loads)
  mapping <- buildMappingObject(
    organism    = organism,
    genomeBuild = genomeBuild,
    txdb        = txdb)


  # mapping <- buildMappingObject("Homo sapiens",getEnsDb = TRUE)

  # 2) Early validate geneType against the OrgDb
  if (!is.null(geneType)) {
    od_cols <- AnnotationDbi::columns(mapping$orgdb)
    if (!"GENETYPE" %in% od_cols) {
      stop(
        "Your OrgDb (", mapping$orgdb, ") does not contain a ‘GENETYPE’ column;\n",
        "cannot apply geneType filter '", geneType, "'.\n",
        "Please omit geneType or choose a supported species OrgDb."
      )
    }
  }

  # 3) Map the user’s IDs — a bit heavier, but now we know geneType is valid
  mapped <- mapForegroundIDs(
    foreground_ids = foreground_ids,
    mapping        = mapping,
    threshold      = threshold,
    transcript     = transcript,
    stripVersions  = stripVersions,
    inflateThresh  = inflateThresh
  )

  # mapped <- mapForegroundIDs(foreground_ids2, mapping, transcript = FALSE)

  # 4) Pool‐level filtering (this is pretty quick, once mapped is in memory)
  filtered <- poolFilter(
    mapped    = mapped,
    geneType  = geneType
  )

  # filtered <- poolFilter(mapped,"protein-coding")


  # # 5) Coordinate extraction (lazy until collect, then quick)
  # coords <- getCoordinates(
  #   mapped      = mapped,
  #   fg_ids      = filtered$fg_ids,
  #   transcript  = transcript,
  #   TSS.method  = ...,
  #   ...
  # )

  # # 6) Finally, the actual enrichment (whatever that is)
  # result <- doEnrichment(coords, …)
  # return(result)

  # Attempt to close the TxDb SQLite handle, if present
  try({
    con <- mapping$so_obj$src$con
    if (!is.null(con)) {
      DBI::dbDisconnect(con)
    }
  }, silent = TRUE)

  return(filtered)

}
#
# transcriptTrialTx <- runDEMotifEnrichment(foreground_ids, transcript = T) # Correctly provides X-row foreground id data table
# transcriptTrialGx <- runDEMotifEnrichment(foreground_ids, transcript = F) # Correctly provides X-row foreground id data table
#
# geneTrialTx <- runDEMotifEnrichment(foreground_ids2, transcript = T) # Correctly errors and asks for different id type
# geneTrialGx <- runDEMotifEnrichment(foreground_ids2, transcript = F)


