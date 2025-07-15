#' @title Build an ID‐mapping table from aspecies-specific OrgDb package
#'
#' @description In response to the user-provided species, install and load the
#' relevant OrgDb package, mapping a base keytype (In most cases ENTREZID) to
#' all other #’ supported identifier columns.
#'
#' @param user_list A vector of identifiers to form the foreground for enrichment
#' analysis
#' @param species string defining species (accept abbreviation or full name)
#' @param threshold Minimal succesful mapping proportion to avoid exiting the function
#' @param verbose boolean for whether informational and progress messages are output
#'
#' @return A tibble with user_ids_, keytype id, mapped (boolean), GENETYPE
#' @export
mapIds <- function(
    user_list,
    stripVersion = TRUE,
    species,
    threshold = 0.9,
    verbose = TRUE) {

    # Clean user_list entries to remove versions if desired
  cleanIds <- function(ids) {
    # pattern: ENSG… .v  OR  ENST… .v OR ENSP… .v  OR  NM_/NR_/NP_ .v
    pat <- "^(ENS[GTPS]\\d+|N[MRP]_\\d+)\\.(\\d+)$"
    idx <- grepl(pat, ids, ignore.case = TRUE)
    ids[idx] <- sub(pat, "\\1", ids[idx], ignore.case = TRUE)
    ids
  }

  if (stripVersion) {
    user_list <- cleanIds(user_list)
  }


  # Identify, install and load the appropriate org.db
  # 1) Normalize the user’s input case
  sp_in <- tolower(species)

  # 2) Test each entry of .speciesMap for a match in code, common or latin
  matches <- vapply(
    .speciesMap,
    FUN = function(info) {
      any(sp_in == tolower(info$code),
          sp_in == tolower(info$common),
          sp_in == tolower(info$latin))
    },
    logical(1)
  )

  # 3) Error if none or more than one hit
  if (!any(matches)) {
    stop("Unknown species ‘", species, "’. Supported are: ",
         paste(names(.speciesMap), collapse = ", "))
  }
  if (sum(matches) > 1) {
    stop("Ambiguous species ‘", species, "’. Try one of: ",
         paste(names(.speciesMap)[matches], collapse = ", "))
  }

  # 4) Pull out the correct entry
  sel_code <- names(.speciesMap)[which(matches)]
  info     <- .speciesMap[[sel_code]]
  orgdb    <- info$orgdb_pkg
  baseKey  <- info$baseKey

  # now install/load orgdb, build map, etc.
  # Report org.db being used
  if (verbose) {
    message(
      "Using OrgDb package '", info$orgdb_pkg,
      "' (", info$latin, ") with base keytype '", info$baseKey, "'."
    )
  }

  # Check if org.db already installed, and if not install it.
  if (!requireNamespace(orgdb, quietly = TRUE)) {
    if (verbose) {
      message("Installing missing OrgDb: ", orgdb)
    }
    BiocManager::install(orgdb, ask = FALSE)
  }

  # Load orgdb package nameSpace so I can access the object
  requireNamespace(orgdb, quietly = TRUE)

  # Get the OrgDb object itself
    #    Most OrgDb pkgs export an object with the same name as the pkg: - need to check if any don't respect this
  od <- get(orgdb, envir = asNamespace(orgdb))

  # list all the columns you want
  all_cols <- AnnotationDbi::columns(od)
  valid_keytypes <- AnnotationDbi::keytypes(od)
  cand_cols      <- intersect(all_cols, valid_keytypes)

  start_time <- Sys.time()
  matched <- vapply(cand_cols, function(col) {
    nmat <- tryCatch({
      df <- AnnotationDbi::select(
        od,
        keys    = user_list,
        keytype = col,
        columns = col
      )
      sum(!is.na(df[[col]]))
    }, error = function(e) {
      # if select() can’t find any keys, just return 0
      0L
    })
    nmat
  }, integer(1))

  end_time <- Sys.time()
  elapsed  <- end_time - start_time
  message("Matching took ", round(as.numeric(elapsed, units="secs"), 2), " seconds")

}

user_list = "SNCA"
mapIds(user_list, species = "hs")

mapIds(user_list, species = "mouse")


