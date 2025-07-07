library(GenomicRanges)
library(GenomicFeatures)  # for makeGRangesFromDataFrame

definePromoterRegions <- function(df, region = c(300, 50)) {
  up   <- abs(region[1])
  down <- region[2]

  # First compute promoter_start/end in the df
  df2 <- df %>%
    mutate(
      promoter_start = if_else(strand == "+", TSS - up,    TSS - down),
      promoter_end   = if_else(strand == "+", TSS + down,  TSS + up)
    )

  # Now make a GRanges, keeping all other columns as metadata
  gr_prom <- makeGRangesFromDataFrame(
    df2,
    keep.extra.columns = TRUE,
    seqnames.field = "seqnames",
    start.field    = "promoter_start",
    end.field      = "promoter_end",
    strand.field   = "strand"
  )

  gr_prom
}
