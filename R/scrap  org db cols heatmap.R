library(AnnotationDbi)
library(dplyr)
library(tidyr)

# assume .speciesMap is already in your namespace
# and has elements of the form list(orgdb_pkg, baseKey, …)

summary_cols <- lapply(names(.speciesMap), function(sp) {
  info <- .speciesMap[[sp]]
  pkg  <- info$orgdb_pkg

  # install if needed, then load namespace
  if (!requireNamespace(pkg, quietly=TRUE))
    BiocManager::install(pkg, ask=FALSE)
  requireNamespace(pkg, quietly=TRUE)

  # extract the OrgDb object
  od <- get(pkg, envir=asNamespace(pkg))

  # get all columns for this OrgDb
  cols <- columns(od)

  # return a tibble row with a list-column of cols
  tibble(
    species    = sp,
    orgdb_pkg  = pkg,
    baseKey    = info$baseKey,
    all_cols   = list(cols)
  )
})

# bind into one tibble
col_df <- bind_rows(summary_cols)

# If you want a long form (one row per species–column pair):
col_long <- col_df %>%
  unnest_longer(all_cols, values_to = "orgdb_column")

# View the wide summary (with list-column)
print(col_df)

# Or view the long version
print(col_long)

library(dplyr)
library(tidyr)
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)   # Bioconductor heatmap package
library(circlize)         # for colorRamp2 if you want more control

mat_df <- col_long %>%
  mutate(present = 1L) %>%             # mark all rows as “present”
  select(species, orgdb_column, present) %>%
  pivot_wider(
    names_from   = column,
    values_from  = present,
    values_fill  = 0L
  )





# 1) Build the long tibble of (species, column)
col_long <- purrr::map_dfr(names(.speciesMap), function(sp) {
  info <- .speciesMap[[sp]]
  pkg  <- info$orgdb_pkg

  if (!requireNamespace(pkg, quietly=TRUE)) {
    BiocManager::install(pkg, ask=FALSE)
  }
  requireNamespace(pkg, quietly=TRUE)

  od <- get(pkg, envir = asNamespace(pkg))
  tibble::tibble(
    species       = sp,
    orgdb_pkg     = pkg,
    baseKey       = info$baseKey,
    orgdb_column  = AnnotationDbi::columns(od)
  )
})

# 2) Pivot to wide 0/1 matrix
mat_df <- col_long %>%
  dplyr::mutate(present = 1L) %>%                            # mark presence
  dplyr::select(species, orgdb_column, present) %>%         # pick the right names
  tidyr::pivot_wider(
    names_from   = orgdb_column,
    values_from  = present,
    values_fill  = list(present = 0L)
  )

# 3) Convert to a matrix (optionally)
mat <- mat_df %>%
  tibble::column_to_rownames("species") %>%
  as.matrix()

# 4) Preview a few rows × cols
knitr::kable(mat[1:6, 1:6],
             caption = "Presence/absence (1/0) of OrgDb columns for first 6 species and 6 columns")

# 5) Draw a simple heatmap with pheatmap
if (requireNamespace("pheatmap", quietly=TRUE)) {
  pheatmap::pheatmap(
    mat,
    cluster_rows    = FALSE,
    cluster_cols    = FALSE,
    legend           = FALSE,
    show_rownames    = TRUE,
    show_colnames    = TRUE
  )
}



install.packages("pheatmap")
library(pheatmap)
# force rendering even inside functions or Rmd chunks

ph <- pheatmap::pheatmap(
  mat,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  legend          = FALSE,
  show_rownames   = TRUE,
  show_colnames   = TRUE
)
print(ph)


colnames(mat)
