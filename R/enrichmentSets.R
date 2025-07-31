#' Run the Full Differential‐Expression Motif Enrichment Pipeline
#'
#' @description
#' \code{enrichmentSets()} is a one‐stop wrapper that takes a vector of
#' user‐provided IDs from a differential expression analysis and returns
#' promoter regions plus matched background sets for downstream motif
#' enrichment.  It performs:
#' 1. ID mapping via \code{\link{buildMappingObject}()} and
#'    \code{\link{mapForegroundIDs}()}.
#' 2. Optional biotype filtering via \code{\link{poolFilter}()}.
#' 3. Promoter coordinate extraction via \code{\link{getCoordinates}()}.
#' 4. Background sampling (pool, random, or matched) within that same call.
#'
#' @param foreground_ids Character vector of user‐supplied gene or transcript IDs
#'   (e.g. Ensembl, RefSeq, gene symbols) to analyze.
#' @param organism Character(1). Species name (e.g. \code{"Homo sapiens"}).
#' @param genomeBuild Character(1). UCSC genome build (e.g. \code{"hg38"}), or
#'   \code{"auto"} to pick the latest supported build.
#' @param txdb Character(1). Name of a \pkg{TxDb} package (e.g.
#'   \code{"TxDb.Hsapiens.UCSC.hg38.knownGene"}), or \code{"auto"}.
#' @param getEnsDb Logical; if \code{TRUE}, also load an \code{EnsDb} for
#'   \code{TSS.method="Ensembl_canonical"}.
#' @param transcript Logical; \code{TRUE} to treat inputs as transcript‐level,
#'   \code{FALSE} for gene‐level.
#' @param threshold Numeric in [0,1]. Min fraction of IDs that must map to
#'   pick a keytype (default 0.9).
#' @param stripVersions Logical; strip “.1”, “.2” suffixes from Ensembl/RefSeq IDs.
#' @param inflateThresh Numeric in [0,1]; max allowed transcript:gene inflation
#'   before auto‐collapsing (default 1).
#' @param geneType Optional character; biotype filter (e.g. \code{"protein-coding"}).
#' @param ensdb Optional \code{EnsDb} object; used only if \code{TSS.method="Ensembl_canonical"}.
#' @param TSS.method Character; TSS selection method for gene‐level mode:
#'   \code{"UCSCgene"}, \code{"Ensembl_canonical"}, \code{"commonTSS"},
#'   \code{"uniqueTSS"}, \code{"fivePrimeTSS"}, or \code{"allTSS"}.
#' @param overlapMinGap Numeric; minimum gap when reducing overlapping promoters.
#' @param onePromoterPerGene Logical; if \code{TRUE}, pick only one promoter per gene.
#' @param bgMethod Character; background sampling method: \code{"matched"},
#'   \code{"pool"}, or \code{"random"}.
#' @param n_ratio Numeric; for \code{bgMethod="random"}, number of bg =
#'   \code{n_ratio * #foreground}.
#' @param bgExcludeFgOverlaps Logical; drop any background overlapping a foreground.
#' @param bgExcludeFgGenes Logical; drop any background gene present in foreground.
#' @param covariates Character vector of covariate names for matching:
#'   \code{"width"}, \code{"gc"}, \code{"chromosome"}, or any custom
#'   \code{mcols()} column.
#' @param bgReplace Logical; allow replacement in sampling (random or matched).
#' @param seed Integer; random seed for reproducibility.
#' @param nrMethod Character; \code{"stratified"}, \code{"rejection"}, or
#'   \code{"nearest"}—passed to matched sampling.
#' @param bsGenome Optional \pkg{BSgenome} object; required if \code{"gc"} is in
#'   \code{covariates} and not otherwise provided.
#' @param promoterWindow Numeric named vector \code{c(upstream, downstream)};
#'   promoter flank widths in bp (default \code{c(300,50)}).
#' @param standardChroms Logical; restrict to standard chromosomes.
#' @param reduceOverlaps Logical; merge overlapping promoter windows.
#'
#' @return A \code{list} containing the output of
#'   \code{\link{getCoordinates}()}, namely:
#'   \describe{
#'     \item{backgroundElements}{\code{GRanges} of sampled background promoters}
#'     \item{foregroundElements}{\code{GRanges} of foreground promoters}
#'     \item{backgroundUniverse}{\code{GRanges} of the pruned promoter pool}
#'     \item{matchObject}{\code{MatchedGRanges} if \code{bgMethod="matched"},
#'       else \code{NULL}}
#'   }
#'
#' @examples
#' \dontrun{
#' # Minimal run with defaults:
#' results <- enrichmentSets(c("ENSG00000139618", "ENSG00000157764"))
#'
#' # Full control:
#' results <- enrichmentSets(
#'   foreground_ids      = my_genes,
#'   organism            = "Homo sapiens",
#'   genomeBuild         = "hg38",
#'   transcript          = TRUE,
#'   geneType            = "protein-coding",
#'   TSS.method          = "commonTSS",
#'   bgMethod            = "matched",
#'   covariates          = c("gc","chromosome"),
#'   nrMethod            = "stratified",
#'   bsGenome            = BSgenome.Hsapiens.UCSC.hg38,
#'   promoterWindow      = c(upstream=500, downstream=100)
#' )
#' }
#'
#' @importFrom AnnotationDbi columns
#' @importFrom Organism.dplyr supportedOrganisms src_organism
#' @importFrom dplyr filter pull
#' @importFrom stringr str_to_lower
#' @importFrom DBI dbDisconnect
#' @export
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



