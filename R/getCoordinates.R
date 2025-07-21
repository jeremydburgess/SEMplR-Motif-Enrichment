getCoordinates <- function(mapped,
                           transcript = FALSE,
                           TSS.method = c(
                             "UCSCgene",
                             "Ensembl_canonical",
                             "commonTSS",
                             "uniqueTSS",
                             "earliestTSS",
                             "allTSS"
                           )) {
  TSS.method <- match.arg(TSS.method)

  # unpack
  fg_df      <- mapped$fg_ids      # data.frame(entrez, mappedID)
  bg_df      <- mapped$bg_ids
  so_obj     <- mapped$so_obj
  transcript <- mapped$transcript   # logical flag
  # (we assume mapped$ensdb already exists, if you need it)
  edb        <- mapped$ensdb

  if (transcript) {
    #
    # === TRANSCRIPT‐LEVEL MODE ===
    #
    # 1) un-collapse mappedIDs → one row per tx_id
    fg_tx <- tidyr::separate_rows(fg_df, mappedID, sep = ",") %>%
      dplyr::rename(tx_id = mappedID)
    bg_tx <- tidyr::separate_rows(bg_df, mappedID, sep = ",") %>%
      dplyr::rename(tx_id = mappedID)

    # 2) pull their full transcript ranges
    tx_ranges <- tbl(so_obj, "ranges_tx") %>%
      dplyr::filter(tx_id %in% bg_tx$tx_id) %>%
      dplyr::select(
        tx_id,
        seqnames = tx_chrom,
        start    = tx_start,
        end      = tx_end,
        strand   = tx_strand
      ) %>%
      collect()

    # 3) restrict to just the foreground transcripts
    fg_ranges <- dplyr::filter(tx_ranges, tx_id %in% fg_tx$tx_id)

    # 4) build GRanges for fg and bg
    gr_bg <- GenomicRanges::makeGRangesFromDataFrame(
      tx_ranges,
      seqnames.field     = "seqnames",
      start.field        = "start",
      end.field          = "end",
      strand.field       = "strand",
      keep.extra.columns = TRUE
    )
    gr_fg <- GenomicRanges::makeGRangesFromDataFrame(
      fg_ranges,
      seqnames.field     = "seqnames",
      start.field        = "start",
      end.field          = "end",
      strand.field       = "strand",
      keep.extra.columns = TRUE
    )

    return(list(bg=gr_bg, fg=gr_fg))
  }

  #
  # === GENE‐LEVEL MODE ===
  #
  if (TSS.method == "UCSCgene") {
    # just join to the gene ranges table
    gene_ranges <- tbl(so_obj, "ranges_gene") %>%
      dplyr::filter(entrez %in% bg_df$entrez) %>%
      dplyr::select(
        entrez,
        seqnames = seqnames,
        start,
        end,
        strand
      ) %>%
      collect()

    fg_ranges <- gene_ranges[gene_ranges$entrez %in% fg_df$entrez, ]

    gr_bg <- GenomicRanges::makeGRangesFromDataFrame(
      gene_ranges, keep.extra.columns=TRUE
    )
    gr_fg <- GenomicRanges::makeGRangesFromDataFrame(
      fg_ranges, keep.extra.columns=TRUE
    )
    return(list(bg=gr_bg, fg=gr_fg))
  }

  # for all other methods (they require transcript‐level TSS computations)
  # 1) un-collapse bg → one tx_id per row, via the id_transcript map
  bg_tx <- tidyr::separate_rows(bg_df, mappedID, sep = ",") %>%
    dplyr::rename(tmpID = mappedID) %>%
    # join to id_transcript to get tx_id
    dplyr::inner_join(
      tbl(so_obj, "id_transcript") %>%
        dplyr::select(entrez, tx_id = ensembltrans) %>%
        collect(),
      by = c("entrez")
    ) %>%
    dplyr::select(entrez, tx_id)

  # 2) if Ensembl_canonical, filter for only canonical transcripts in EnsDb
  if (TSS.method == "Ensembl_canonical") {
    can_tx <- ensembldb::transcripts(edb) %>%
      as.data.frame() %>%
      dplyr::filter(tx_is_canonical == TRUE) %>%
      dplyr::pull(tx_id)
    bg_tx <- dplyr::filter(bg_tx, tx_id %in% can_tx)
  }

  # 3) fetch full ranges for the (possibly filtered) tx_ids
  tx_coords <- tbl(so_obj, "ranges_tx") %>%
    dplyr::filter(tx_id %in% bg_tx$tx_id) %>%
    dplyr::select(
      tx_id,
      seqnames = tx_chrom,
      start    = tx_start,
      end      = tx_end,
      strand   = tx_strand
    ) %>%
    collect() %>%
    # re‐attach entrez
    dplyr::inner_join(bg_tx, by = "tx_id") %>%
    dplyr::mutate(tss = ifelse(strand == "+", start, end))

  # 4) apply the TSS‐method collapse (skip for 'allTSS')
  tss_df <- switch(
    TSS.method,
    commonTSS   = tx_coords %>% dplyr::group_by(entrez) %>%
      dplyr::slice_max(tabulate(tss), with_ties=TRUE), # keeps multiple if more than one max TSS frequency
    uniqueTSS   = tx_coords %>% dplyr::group_by(entrez,tss) %>% dplyr::slice(1),
    earliestTSS = tx_coords %>% dplyr::group_by(entrez) %>%
      dplyr::slice_min(tss, with_ties=FALSE) # Only keeps,one of tied results
    allTSS      = tx_coords
  )

  # 5) build GRanges
  gr_bg <- GenomicRanges::GRanges(
    seqnames = tss_df$seqnames,
    ranges   = IRanges::IRanges(start=tss_df$start
                                end=tss_df$end),
    strand   = tss_df$strand,
    tx_id    = tss_df$tx_id,
    entrez   = tss_df$entrez,
    tss      = tss_df$tss
  )
  # filter FG down to those in fg_df
  fg_ids_split <- tidyr::separate_rows(fg_df, mappedID, sep=",") %>%
    dplyr::rename(tx_id = mappedID)
  gr_fg <- gr_bg[gr_bg$tx_id %in% fg_ids_split$tx_id]

  return(list(bg=gr_bg, fg=gr_fg))
}
