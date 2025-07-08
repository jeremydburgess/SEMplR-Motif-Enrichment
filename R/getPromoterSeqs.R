library(GenomicRanges)
library(Biostrings)

getPromoterSeqs <- function(
    prom_gr,          # GRanges, or named list of GRanges, or GRangesList
    bsgenome,         # e.g. BSgenome.Hsapiens.UCSC.hg38
    flank_rc = TRUE   # whether to reverse‐complement minus‐strand sequences
) {
  # 1) If it's a GRangesList or list, recurse over each element
  if (is(prom_gr, "GRangesList") || is.list(prom_gr)) {
    out <- lapply(prom_gr, function(gr) {
      getPromoterSeqs(gr, bsgenome, flank_rc)
    })
    names(out) <- names(prom_gr)

    # If input was a GRangesList, wrap into a DNAStringSetList
    if (is(prom_gr, "GRangesList")) {
      out <- do.call(DNAStringSetList, out)
      names(out) <- names(prom_gr)
    }
    return(out)
  }

  # 2) Otherwise prom_gr is a single GRanges → getSeq + RC on “-” strands
  seqs <- getSeq(bsgenome, prom_gr)
  if (flank_rc) {
    is_neg <- as.vector(strand(prom_gr) == "-")
    seqs[is_neg] <- reverseComplement(seqs[is_neg])
  }

  # 3) Name by your unique identifier column
  names(seqs) <- mcols(prom_gr)$unique_id

  seqs
}
