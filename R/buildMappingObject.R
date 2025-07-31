#' Build an organism-specific id-mapping object
#'
#' @description
#' \code{buildMappingObject()} installs (if needed) and loads the appropriate
#' \pkg{OrgDb} and \pkg{TxDb} Bioconductor packages for a given species and
#' genome build, then returns an \code{src_organism} object along with metadata
#' about which packages and parameters were used.
#'
#' @param organism Character scalar specifying the species, e.g. \dQuote{Homo sapiens}.
#' @param genomeBuild Character scalar giving the genome build (e.g. \dQuote{hg38}),
#'   or \dQuote{auto} (default) to pick the latest supported build.
#' @param txdb Character scalar naming a specific \pkg{TxDb} package, or
#'   \dQuote{auto} (default) to select the most recent \pkg{TxDb} for the build.
#'
#' @return A named list with components:
#' \describe{
#'   \item{so_obj}{An \code{src_organism} data source (from \pkg{Organism.dplyr}).}
#'   \item{orgdb}{The loaded \pkg{OrgDb} package namespace object.}
#'   \item{organism}{The species string used.}
#'   \item{genomeBuild}{The resolved genome build string.}
#'   \item{txdb}{The resolved \pkg{TxDb} package name.}
#' }
#'
#' @examples
#' \dontrun{
#' mapping <- buildMappingObject(
#'   organism    = "Homo sapiens",
#'   genomeBuild = "auto",
#'   txdb        = "auto"
#' )
#' so_obj <- mapping$so_obj
#' }
#'
#' @seealso \code{\link[Organism.dplyr]{supportedOrganisms}},
#'   \code{\link[Organism.dplyr]{src_organism}}
#'
#' @importFrom Organism.dplyr supportedOrganisms src_organism
#' @importFrom stringr str_to_lower
#' @importFrom dplyr filter pull
#' @importFrom jamba mixedSort
#' @importFrom magrittr %>%
#' @export
buildMappingObject <- function(
    organism    = "Homo sapiens",
    genomeBuild = "auto",       # “auto” → pick latest for that organism
    txdb        = "auto",      # “auto” → pick latest for that build
    getEnsDb    = FALSE)      # Also load EnDb (used for Ensembl_canonical filtering later)
 {

  # 1) figure out which TxDb’s are supported for this organism
  so <- Organism.dplyr::supportedOrganisms()
  hits <- so %>% dplyr::filter(str_to_lower(organism) == str_to_lower(!!organism))
  if (nrow(hits)==0) {
    stop("Organism ‘", organism, "’ is not supported. Supported organisms are: ",
         paste0(unique(so$organism),collapse = ", "))
  }

  # Store required OrgDb
  orgdb <- hits$OrgDb[1]
  message("Using orgDb: ",orgdb)

  # 2) Extract the set of available genomeBuilds for that organism
  builds <- unique(hits$TxDb %>%
                     sub("^TxDb\\.[^.]+\\.UCSC\\.", "", .) %>%
                     sub("\\..*$", "", .))

  # 3) resolve genomeBuild
  if (is.null(genomeBuild) || genomeBuild=="auto") {
    # pick the “latest” build by lexical sort (or you could parse numbers)
    genomeBuild <- jamba::mixedSort(builds)[length(builds)]
  } else if (!(genomeBuild %in% builds)) {
    stop(
      "Invalid genomeBuild '", genomeBuild,
      "' for organism '", organism, "'.\n",
      "Available builds are: ", paste(builds, collapse=", ")
    )
  }

  message("Using GenomeBuild: ",genomeBuild)

  # 4) narrow down to the matching TxDb records
  tx_choices <- hits %>%
    dplyr::filter(grepl(paste0("\\.", genomeBuild, "\\."), TxDb)) %>%
    pull(TxDb)
  if (length(tx_choices)==0) {
    stop("No TxDb found for ", organism, " + build=", genomeBuild)
  }

  # 5) pick the ‘txdb’ if user asked for a specific one
  if (txdb!="auto") {
    if (!(txdb %in% tx_choices)) {
      stop("Invalid txdb '", txdb, "'. Available for this build:\n",
           paste(tx_choices, collapse="\n"))
    }
    tx_choice <- txdb
  } else {
    tx_choice <- sort(tx_choices)[length(tx_choices)]
  }

  message("using TxDb: ",tx_choice)

  # If necessary, check for and load EnsDb
  ensdb <- NULL
  if(isTRUE(getEnsDb)) {
    ensdb <- helper_loadEnsDbForOrganism(organism,genomeBuild)}

  # Check both orgDb and TxDbs are installed and call them into the memory
  # ensure packages are loaded
  # later in your wrapper, after you’ve resolved orgdb and tx_choice:
  helper_ensureInstalledAndLoadedBioC(orgdb)
  orgdb_obj <- get(orgdb, envir = asNamespace(orgdb))

  helper_ensureInstalledAndLoadedBioC(tx_choice)


  # 6) now you can do
  so_obj <- Organism.dplyr::src_organism(
    txdb     = tx_choice
  )

  so_obj <- Organism.dplyr::src_organism(txdb = tx_choice)

  # Return *both* the data‐source objects *and* the parameters that created them
  list(
    so_obj      = so_obj,
    orgdb       = orgdb_obj,    # the actual OrgDb, not its name
    organism    = organism,
    genomeBuild = genomeBuild,
    txdb        = tx_choice,
    ensdb       = ensdb
  )

  }
