#' Reduce overlapping ranges within each gene
#'
#' Performs a two-step reduction of genomic ranges grouped by gene:
#'
#' 1. A global reduction across all ranges using `extraChIPs::reduceMC()` to quickly identify
#'    any promoter windows that merge across different genes.
#' 2. A per-gene reduction for only those genes whose windows were merged in the global pass,
#'    by subsetting the original ranges, splitting by gene ID, and running `reduceMC()` on each subset.
#'
#' @param gr A `GRanges` object containing promoter ranges with an `entrez` metadata column.
#' @param overlapMinGap Integer(1); the `min.gapwidth` argument passed to `reduceMC()` (default 0).
#' @return A `GRanges` object of reduced promoter windows, where overlaps are merged only within each gene,
#' preserving all metadata columns.
#' @import GenomicRanges
#' @importFrom extraChIPs reduceMC
#' @importFrom BiocGenerics sort
#' @importFrom pbapply pblapply
#' @export
helper_reduceOverlapsWithinGenes <- function(gr, overlapMinGap = 0) {
  # 1) Global reduction (may merge across genes)
  red_global <- extraChIPs::reduceMC(gr, min.gapwidth = overlapMinGap)

  # 2) Identify genes whose windows merged across genes
  entrez_list <- mcols(red_global)$entrez
  mixed_genes <- unique(unlist(entrez_list[lengths(entrez_list) > 1]))
  all_genes   <- unique(mcols(gr)$entrez)
  safe_genes  <- setdiff(all_genes, mixed_genes)

  # 3) Subset the original input by safe vs mixed genes
  orig_safe  <- gr[mcols(gr)$entrez %in% safe_genes]
  orig_mixed <- gr[mcols(gr)$entrez %in% mixed_genes]

  # 4) Single global reduce for safe genes
  red_safe <- extraChIPs::reduceMC(orig_safe,
                                   min.gapwidth = overlapMinGap)

  # 5) Split and reduce per gene for mixed genes
  grp_mixed <- split(orig_mixed, mcols(orig_mixed)$entrez)
  red_mixed_list <- pbapply::pblapply(grp_mixed, function(subgr) {
    extraChIPs::reduceMC(subgr, min.gapwidth = overlapMinGap)
  })

  # Combine the per-gene reductions into one GRanges
  red_mixed <- GenomicRanges::GRangesList(red_mixed_list) %>%
    unlist(use.names = FALSE)

  # 6) Combine safe and mixed results, then sort
  final_promoters <- BiocGenerics::sort(c(red_safe, red_mixed))

  return(final_promoters)
}

