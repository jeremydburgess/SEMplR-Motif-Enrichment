#' Filter Foreground and Background ID Sets by Gene Type
#'
#' @description
#' `poolFilter()` takes the mapped foreground and background ID data frames
#' (as produced by `mapIDs()`) and, if requested, filters both sets
#' to only include genes (or their transcripts) of a specified biotype
#' (e.g. “protein-coding”).
#'
#' @param mapped    A list returned by `mapIDs()`, containing at least:
#'   \itemize{
#'     \item `fg_ids`: data.frame with columns `entrez` and `mappedID` (foreground).
#'     \item `bg_ids`: data.frame with columns `entrez` and `mappedID` (background).
#'     \item `so_obj`: a `src_organism` object for transcript lookups.
#'     \item `orgdb`:  the loaded OrgDb package object.
#'     \item `transcript`: logical, whether IDs are transcripts.
#'   }
#' @param geneType  Optional character scalar; if not NULL, only genes of this
#'   biotype (`GENETYPE` in the OrgDb) will be retained.  Valid values vary by
#'   organism (e.g. “protein-coding”, “lncRNA”, etc.).
#'
#' @return
#' The original `mapped` list, but with `fg_ids` and `bg_ids` replaced by
#' filtered versions (only rows matching `geneType`, if provided).
#'
#' @examples
#' \dontrun{
#' mapping  <- buildMappingObject("Homo sapiens")
#' mapped   <- mapIDs(my_ids, mapping)
#' filtered <- poolFilter(mapped, geneType = "protein-coding")
#' }
#'
#' @importFrom dplyr tbl filter select collect
#' @keywords internal
poolFilter <- function(mapped,
                       geneType = NULL) {


  # unpack
  fg_df      <- mapped$fg_ids
  bg_df      <- mapped$bg_ids
  transcript <- mapped$transcript
  so_obj     <- mapped$so_obj
  orgdb      <- mapped$orgdb

  # geneType filtering, only if requested
  if (!is.null(geneType)) {
    # fetch every gene → GENETYPE mapping for the background universe
    gt <- suppressMessages(
      AnnotationDbi::select(
        orgdb,
        keys    = bg_df$entrez,
        keytype = "ENTREZID",
        columns = c("ENTREZID", "GENETYPE")
      )
    )
    valid_types <- sort(unique(gt$GENETYPE))

    # if the user’s geneType isn’t in that set, stop and list the valid ones
    if (!(geneType %in% valid_types)) {
      stop(
        "Invalid geneType '", geneType, "'.\n",
        "Valid geneType values in this OrgDb are:\n  ",
        paste(valid_types, collapse = ", ")
      )
    }

    # now subset to the requested geneType
    keep_genes <- gt$ENTREZID[gt$GENETYPE == geneType]

    if (transcript) {
      # transcripts → genes map
      idtx <- tbl(so_obj, "id_transcript") %>%
        dplyr::filter(entrez %in% keep_genes) %>%
        dplyr::select(entrez, ensembltrans) %>%
        collect()

      # restrict both bg and fg to only those transcripts whose gene is in keep_genes
      bg_df <- bg_df[
        bg_df$entrez   %in% keep_genes &
          bg_df$mappedID %in% idtx$ensembltrans,
      ]
      fg_df <- fg_df[
        fg_df$mappedID %in% idtx$ensembltrans,
      ]

    } else {
      # gene‐mode: intersect by entrez
      bg_df <- bg_df[ bg_df$entrez %in% keep_genes, ]
      fg_df <- fg_df[ fg_df$entrez %in% keep_genes, ]
    }
  }

  # reinject filtered dfs and return all components
  mapped$fg_ids <- fg_df
  mapped$bg_ids <- bg_df

  return(mapped)
}


