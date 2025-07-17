# R/choosePath.R
# Function to determine the available arguments and path taken for mapping, filtering,
# and enriching differentially expressed genes

choosePath() <- function(species = "Hs",
                         genomeBuild = "auto",
                         ensemblRelease = "auto",
                         verbose = TRUE) {

  # Resolve species → latin name for matching to EnsDbs in Annotation hub
  sp_in <- tolower(species)
  hits  <- vapply(.speciesMap, function(info) {
    any(sp_in == tolower(info$code),
        sp_in == tolower(info$common),
        sp_in == tolower(info$latin))
  }, logical(1))
  if (!any(hits))  stop("Unknown species: ", species)
  if (sum(hits)>1) stop("Ambiguous species: ", species)

  latinSpecies <- info$latin

  # Resolve latin name → format first letter of GenusSpecies (eg 'Hsapiens') for
  # matching to TxDbs
  latin_to_txdb <- function(latin_name) {
    parts <- strsplit(latin_name, "\\s+")[[1]]
    if (length(parts) < 2) {
      stop("'", latin_name, "' is not a two-part Latin name")
    }
    genus   <- parts[1]
    species <- parts[2]
    paste0(toupper(substr(genus, 1, 1)), tolower(species))
  }

  TxDbSpecies <- latin_to_txdb(latinSpecies)




}

BiocManager::install(c("Organism.dplyr"))
library(Organism.dplyr)
# 1. See what organisms & TxDbs are supported out of the box
supported <- supportedOrganisms()
print(supported, n = Inf)

unique(supported$organism)


unique(ensd_113$species)


supportedOrganisms()
