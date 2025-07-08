library(GenomicRanges)
library(GenomicFeatures)  # for makeGRangesFromDataFrame
library(GenomicRanges)

library(dplyr)
library(GenomicRanges)
library(GenomicFeatures)

definePromoterRegions <- function(pools,
                                  region = c(upstream = 300, downstream = 50)) {
  up   <- abs(region["upstream"])
  down <- region["downstream"]

  GRangesList(
    backgroundUniverse = {
      df2 <- pools$backgroundUniverse %>%
        mutate(
          promoter_start = if_else(strand == "+", TSS - up,    TSS - down),
          promoter_end   = if_else(strand == "+", TSS + down,  TSS + up)
        )
      makeGRangesFromDataFrame(df2, keep.extra.columns=TRUE,
                               seqnames.field="seqnames", start.field="promoter_start",
                               end.field="promoter_end", strand.field="strand")
    },
    foregroundElements = {
      df2 <- pools$foregroundElements %>%
        mutate(
          promoter_start = if_else(strand == "+", TSS - up,    TSS - down),
          promoter_end   = if_else(strand == "+", TSS + down,  TSS + up)
        )
      makeGRangesFromDataFrame(df2, keep.extra.columns=TRUE,
                               seqnames.field="seqnames", start.field="promoter_start",
                               end.field="promoter_end", strand.field="strand")
    }
  )
}

# usage:
# prom_grl <- definePromoterRegions(pools)
# prom_grl$backgroundUniverse  # GRanges
# prom_grl$foregroundElements  # GRanges
