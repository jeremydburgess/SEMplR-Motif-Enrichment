getCoordinates <- function(mapped,
                           fg_ids_df,
                           bg_ids_df,
                           method = c(
                             "UCSCGene",
                             "EnsemblCanonical",
                             "MANE_Select",
                             "commonTSS",
                             "uniqueTSS",
                             "allTSS",
                             "earliestTSS"
                           ),
                           ensdb  = NULL
) {
  method     <- match.arg(method)
  so_obj     <- mapped$so_obj
  transcript <- mapped$transcript

  # Internal helper: take a vector of IDs and produce GRanges
  coords_for <- function(ids) {
    # Gene vs transcript branch
    if (!transcript) {
      # UCSCGene is just ranges_gene
      if (method == "UCSCGene") {
        df <- tbl(so_obj, "ranges_gene") %>%
          filter(id %in% ids) %>%
          select(id, seqnames, start, end, strand) %>%
          collect()
        return(GenomicRanges::makeGRangesFromDataFrame(
          df,
          seqnames.field     = "seqnames",
          start.field        = "start",
          end.field          = "end",
          strand.field       = "strand",
          keep.extra.columns = TRUE
        ))
      }
      # All other methods require transcripts → TSS logic
      tx_ids <- ids
    } else {
      tx_ids <- ids
    }

    # If we get here, method is either EnsemblCanonical/MANE (which handle genes via ensdb)
    # or one of the TSS methods.

    if (method %in% c("EnsemblCanonical","MANE_Select")) {
      if (is.null(ensdb)) {
        stop("For method '", method, "' you must supply an `ensdb` object.")
      }
      txs <- ensembldb::transcripts(
        ensdb,
        columns = c("tx_id","gene_id","mane_select"),
        filter  = TxIdFilter(tx_ids)
      )
      if (method == "EnsemblCanonical") {
        sel <- ensembldb::selectCanonical(txs)
      } else {
        sel <- txs[txs$mane_select == TRUE, ]
      }
      return(GenomicRanges::granges(sel, use.mcols=TRUE))
    }

    # Otherwise we’re in one of the TSS methods
    tx_df <- tbl(so_obj, "ranges_tx") %>%
      filter(tx_id %in% tx_ids) %>%
      select(
        tx_id,
        seqnames = tx_chrom,
        start    = tx_start,
        end      = tx_end,
        strand   = tx_strand
      ) %>%
      collect() %>%
      mutate(tss = if_else(strand == "+", start, end))

    if (method != "allTSS") {
      idtx <- tbl(so_obj, "id_transcript") %>%
        filter(ensembltrans %in% tx_df$tx_id) %>%
        select(entrez, ensembltrans) %>%
        collect()
      tx_df <- left_join(tx_df, idtx, by = c("tx_id" = "ensembltrans"))
    }

    tss_df <- switch(
      method,
      allTSS      = tx_df,
      commonTSS   = {
        counts <- tx_df %>% count(entrez, tss)
        most   <- counts %>% group_by(entrez) %>% slice_max(n, with_ties=FALSE)
        inner_join(most, tx_df, by=c("entrez","tss"))
      },
      uniqueTSS   = tx_df %>% distinct(entrez,tss,.keep_all=TRUE),
      earliestTSS = tx_df %>% group_by(entrez) %>% slice_min(tss, with_ties=FALSE)
    )

    GenomicRanges::GRanges(
      seqnames = tss_df$seqnames,
      ranges   = IRanges::IRanges(start = tss_df$tss, width = 1),
      strand   = tss_df$strand,
      tx_id    = tss_df$tx_id,
      entrez   = tss_df$entrez,
      tss      = tss_df$tss
    )
  }

  # now apply to both sets
  fg_gr <- coords_for(if (!transcript) fg_ids_df$entrez else fg_ids_df$mappedID)
  bg_gr <- coords_for(if (!transcript) bg_ids_df$entrez else bg_ids_df$mappedID)

  # return both
  list(
    fg_gr = fg_gr,
    bg_gr = bg_gr
  )
}
