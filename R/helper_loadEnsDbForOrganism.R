#' Find and load the appropriate EnsDb for a given organism + genome build
#'
#' @description
#' \code{helper_loadEnsDbForOrganism()} looks up the correct Ensembl
#' transcript database (EnsDb) in AnnotationHub for a given species
#' and UCSC‐style genome build (e.g. “hg38”), and returns the freshest
#' EnsDb object.
#'
#' @param organism    Character(1). Latin or common name of the organism
#'                    (e.g. “Homo sapiens” or “human”).
#' @param genomeBuild Character(1). UCSC‐style genome build (e.g. “hg38”).
#'
#' @return An \code{EnsDb} object loaded from AnnotationHub.
#'
#' @details
#' 1. Uses \code{GenomeInfoDb::mapGenomeBuilds()} to translate the UCSC
#'    build into one or more Ensembl build identifiers (including patched
#'    versions).
#' 2. Queries \code{AnnotationHub()} for all “EnsDb” records for the
#'    requested organism.
#' 3. Filters those in‐memory to retain only the builds matching your
#'    UCSC → Ensembl mapping.
#' 4. Picks the record with the latest \code{rdatadateadded} metadata and
#'    returns it.
#'
#' @importFrom GenomeInfoDb mapGenomeBuilds
#' @importFrom AnnotationHub AnnotationHub query subset
#' @keywords internal
helper_loadEnsDbForOrganism <- function(organism, genomeBuild) {
  # 1) Fetch UCSC→Ensembl build mapping
  #    mapGenomeBuilds lives in GenomeInfoDb (Bioc), returns both IDs.
  gbmap <- GenomeInfoDb::mapGenomeBuilds(genomeBuild, style="Ensembl")
  if (nrow(gbmap) == 0) {
    stop("Unknown genomeBuild ‘", genomeBuild, "’")
  }
  # Ensembl genome builds are patched so have 1 UCSC -> multiple ensembl style genome builds.
  # We need to pass the full list by collapsing into OR-regex
  buildPattern <- paste(gbmap$ensemblID, collapse="|")

  # 2) Query AnnotationHub for EnsDb records matching organism + build
  ah <- AnnotationHub()
  speciesHits <- query(ah, c("EnsDb", organism))

  if (length(speciesHits) == 0) {
    stop("No EnsDb available for ", organism)
  }

  # 3) Restrict by genome
  hits <- query(speciesHits,buildPattern)

  if (length(hits) == 0) {
    stop("No EnsDb available for ", genomeBuild)
  }

  # 3) pick the most recent by rdatadateadded
  latestMatch <- AnnotationHub::subset(hits,rdatadateadded == max(hits$rdatadateadded))[1]

  # 4) load & return the EnsDb
  ensdb <- ah[[latestMatch$ah_id]]
  return(ensdb)
}
