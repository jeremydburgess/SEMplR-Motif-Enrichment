poolFilter <- function(mapped,
                       geneType = NULL) {
  # mapped: list from mapForegroundIDs(), containing
  #   fg_ids     — data.frame(entrez, mappedID)
  #   bg_ids     — data.frame(entrez, mappedID)
  #   userIDtype — string
  #   transcript — logical
  #   so_obj     — src_organism
  #   orgdb      — OrgDb package

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
        filter(entrez %in% keep_genes) %>%
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
