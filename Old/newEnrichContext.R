# Constructor for pipeline context
newEnrichContext <- function(species, genomeBuild = "auto") {
  ctx <- list(
    species     = species,
    genomeBuild = genomeBuild,
    soObj       = NULL,      # placeholder for src_organism
    fgEnsembl   = NULL,      # placeholder for mapped foreground IDs
    results     = NULL       # placeholder for final enrichment output
  )
  class(ctx) <- "EnrichContext"
  ctx
}
