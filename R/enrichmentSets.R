enrichmentSets <- function(foreground_ids, # mapForegroundIds
                           organism              = "Homo sapiens", # buildMappingObject, poolFilter, getCoordinates
                           genomeBuild           = "auto", # buildMappingObject, poolFilter, getCoordinates
                           txdb                  = "auto", # buildMappingObject
                           getEnsDb              = FALSE, # buildMappingObject
                           transcript            = FALSE, # mapForegroundIDs, getCoordinates
                           threshold             = 0.9, # mapForegroundIDs
                           stripVersions         = TRUE, # mapForegroundIDs
                           inflateThresh         = 1, # mapForegroundIDs
                           geneType              = NULL, # poolFilter
                           ensdb                 = NULL, # poolFilter
                           TSS.method            = c("UCSCgene","Ensembl_canonical","commonTSS","uniqueTSS","fivePrimeTSS","allTSS"), # getCoordinates
                           overlapMinGap         = 0, # getCoordinates
                           onePromoterPerGene    = FALSE, # getCoordinates
                           bgMethod              = c("matched","pool","random"), # getCoordinates
                           n_ratio               = 1, # getCoordinates
                           bgExcludeFgOverlaps   = TRUE, # getCoordinates
                           bgExcludeFgGenes      = FALSE, # getCoordinates
                           covariates            = c("width", "gc", "chromosome"), # getCoordinates
                           bgReplace             = FALSE, # getCoordinates
                           seed                  = 1985, # getCoordinates
                           nrMethod              = c("stratified","rejection","nearest"), # getCoordinates
                           bsGenome              = NULL, # getCoordinates
                           promoterWindow        = c(upstream=300, downstream=50), # getCoordinates
                           standardChroms        = TRUE, # getCoordinates
                           reduceOverlaps        = TRUE  # getCoordinates
                            ) {

  TSS.method <- match.arg(TSS.method)
  bgMethod <- match.arg(bgMethod)
  nrMethod <- match.arg(nrMethod)

  # early loading of ensdb if required
  if(TSS.method == "Ensembl_canonical" && is.null(ensdb)){getEnsDb <- TRUE}

  # 1) Build the mapping object (cheap metadata + package loads)
  mapping <- buildMappingObject(
    organism    = organism,
    genomeBuild = genomeBuild,
    txdb        = txdb,
    getEnsDb = getEnsDb)

  # For testing:
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

  # For testing:
  # mapped <- mapForegroundIDs(foreground_ids2, mapping, transcript = FALSE)

  # 4) Pool‐level filtering (this is pretty quick, once mapped is in memory)
  filtered <- poolFilter(
    mapped    = mapped,
    geneType  = geneType
  )

  # For testing:
  # filtered <- poolFilter(mapped,"protein-coding")

  # # 5) Coordinate extraction (lazy until collect, then quick)
  coords <- getCoordinates(
    mapped = filtered,
    transcript = transcript,
    TSS.method = TSS.method,
    bgMethod             = bgMethod,
    n_ratio              = n_ratio,
    bgExcludeFgOverlaps  = bgExcludeFgOverlaps,
    bgExcludeFgGenes     = bgExcludeFgGenes,
    covariates           = covariates,
    bgReplace            = bgReplace,
    seed                 = seed,
    nrMethod             = nrMethod,
    bsGenome             = bsGenome,
    promoterWindow       = promoterWindow,
    standardChroms       = standardChroms,
    reduceOverlaps       = reduceOverlaps,
    overlapMinGap        = overlapMinGap,
    onePromoterPerGene   = onePromoterPerGene,
    ensdb                = ensdb)

  # Attempt to close the TxDb SQLite handle, if present
  try({
    con <- mapping$so_obj$src$con
    if (!is.null(con)) {
      DBI::dbDisconnect(con)
    }
  }, silent = TRUE)

  return(coords)

}
#
# transcriptTrialTx <- runDEMotifEnrichment(foreground_ids, transcript = T) # Correctly provides X-row foreground id data table
# transcriptTrialGx <- runDEMotifEnrichment(foreground_ids, transcript = F) # Correctly provides X-row foreground id data table
#
# geneTrialTx <- runDEMotifEnrichment(foreground_ids2, transcript = T) # Correctly errors and asks for different id type
# geneTrialGx <- runDEMotifEnrichment(foreground_ids2, transcript = F)


