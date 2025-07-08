extractPromoterSeqs <- function(
    prom_gr,          # GRanges of promoter windows
    bsgenome,         # e.g. BSgenome.Hsapiens.UCSC.hg38
    flank_rc = TRUE   # whether to reverse‐complement minus‐strand sequences
) {
  seqs <- getSeq(bsgenome, prom_gr)
  if (flank_rc) {
    # reverse‐complement any ranges on the negative strand
    is_neg <- as.vector(strand(prom_gr) == "-")
    seqs[is_neg] <- reverseComplement(seqs[is_neg])
  }
  # attach back as a DNAStringSet with names linking to metadata
  names(seqs) <- mcols(prom_gr)$unique_id
  seqs
}
