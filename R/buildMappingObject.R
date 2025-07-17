# install.packages("jamba")
# library(jamba)

buildMappingObject <- function(
    organism = "Homo sapiens",
    genomeBuild = "auto",       # or NULL
    txdb         = "auto")      # “auto” → pick latest for that build
 {

              # Function to install required bioconductor packages (orgDbs and txDbs)
              ensure_installed_and_loaded <- function(pkg) {
              if (!requireNamespace(pkg, quietly = TRUE)) {
                message("Installing Bioconductor package ‘", pkg, "’ …")
                BiocManager::install(pkg, ask = FALSE, update = FALSE)
              }
              # now load it
              library(pkg, character.only = TRUE)
            }

    # 1) figure out which TxDb’s are supported for this species
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

  # Check both orgDb and TxDbs are installed and call them into the memory
  # ensure packages are loaded
  # later in your wrapper, after you’ve resolved orgdb and tx_choice:
  ensure_installed_and_loaded(orgdb)
  ensure_installed_and_loaded(tx_choice)


  # 6) now you can do
  so_obj <- Organism.dplyr::src_organism(
    txdb     = tx_choice
  )

  return(so_obj)
}
