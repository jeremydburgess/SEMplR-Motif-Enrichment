enrichDEGMotifs <- function(foreground_ids, # vector of ids of genes/transcripts to be used as enrichment foreground
                            species, # species of source organism
                            genomeBuild = "auto",
                            ){

  ctx <- newEnrichContext(species = species,
                          genomeBuild = genomeBuild)





}

so_obj <- buildMappingObject(organism = "Homo sapiens")
mapForegroundIDs(foreground_ids = foreground_ids,so_obj = so_obj)



library(stringr)
library(jamba)
library(tidyr)
library(dplyr)




library(dplyr)

library(dplyr)

# 1. Pull all unique Ensembl gene IDs
all_ensembltrans_ids <- tbl(so_obj, "id_transcript") %>%
  dplyr::select(ensembltrans) %>%
  distinct() %>%
  collect() %>%
  pull(ensembltrans)

# 2. Sample up to 500 of them
set.seed(123)
foreground_ids <- sample(
  all_ensembltrans_ids,
  size    = min(500, length(all_ensembl_ids)),
  replace = FALSE
)

# inspect
length(foreground_ids)  # should be 500 (or fewer if the table was smaller)
head(foreground_ids)


src_tbls(so_obj)
genes_tbl <- tbl(so_obj, "ranges_gene")



safe_map_count <- function(src, keys, kt) {
  if (!(kt %in% valid_kts)) return(NA_integer_)
  result <- tryCatch(
    AnnotationDbi::mapIds(
      src,
      keys     = keys,
      column   = kt,
      keytype  = kt,
      multiVals = "first"
    ),
    error = function(e) rep(NA_character_, length(keys))
  )
  sum(!is.na(result))
}


stats_fast <- data.frame(
  id_type = fast_cols,
  matched = vapply(
    fast_cols,
    function(kt) safe_map_count(so_obj, foreground_ids, kt),
    integer(1)
  ),
  stringsAsFactors = FALSE
)
stats_fast$pct_matched <- stats_fast$matched / length(foreground_ids) * 100
