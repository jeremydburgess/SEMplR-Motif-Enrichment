library(nullranges)    # for matchRanges()
library(Biostrings)    # for getSeq(), letterFrequency()
library(GenomicRanges)

defineBackgroundElements <- function(
    background_universe,     # GRanges
    foreground_elements,     # GRanges
    bsgenome,                # BSgenome object, e.g. BSgenome.Hsapiens.UCSC.hg38
    method            = c("pool", "random", "matched"),
    n_ratio           = 1,
    covariates        = c("width", "gc"),
    exclude_fg        = TRUE,
    allow_replacement = FALSE,
    seed              = NULL,
    nullranges_args   = list(method = "rejection", replace = FALSE, verbose = FALSE)
) {
  method <- match.arg(method)

  # (1) Optionally drop the true foreground from the universe
  if (exclude_fg) {
    background_universe <- subsetByOverlaps(
      background_universe, foreground_elements, invert = TRUE
    )
  }

  # (2) “pool” or “random” are trivial
  if (method == "pool") {
    return(background_universe)
  }
  if (method == "random") {
    if (!is.null(seed)) set.seed(seed)
    n_bg <- n_ratio * length(foreground_elements)
    idx <- sample(
      seq_along(background_universe),
      size    = n_bg,
      replace = allow_replacement
    )
    return(background_universe[idx])
  }

  # (3) matched sampling via nullranges::matchRanges()
  #    — ensure covariate columns exist in mcols(background_universe) & foreground
  #    — width and gc are common choices

  # 3a) width
  if ("width" %in% covariates) {
    # only add if not already present
    if (!"width" %in% names(mcols(background_universe))) {
      mcols(background_universe)$width    <- width(background_universe)
      mcols(foreground_elements)$width    <- width(foreground_elements)
    }
  }

  # 3b) GC content
  if ("gc" %in% covariates) {
    # compute GC fraction per range
    seqs_bg <- getSeq(bsgenome, background_universe)
    seqs_fg <- getSeq(bsgenome, foreground_elements)
    mcols(background_universe)$gc    <- letterFrequency(seqs_bg, letters="GC", as.prob=TRUE)
    mcols(foreground_elements)$gc    <- letterFrequency(seqs_fg, letters="GC", as.prob=TRUE)
  }

  # 3c) build the formula for matchRanges()
  covar_formula <- as.formula(
    paste0("~", paste(covariates, collapse = " + "))
  )

  # set character covariates to factors
  for (cov in covariates) {
    col_data_bg <- mcols(background_universe)[[cov]]
    col_data_fg <- mcols(foreground_elements)[[cov]]

    # if it’s a character vector, make it a factor (preserving unique levels)
    if (is.character(col_data_bg)) {
      mcols(background_universe)[[cov]] <- factor(col_data_bg)
    }
    if (is.character(col_data_fg)) {
      mcols(foreground_elements)[[cov]] <- factor(col_data_fg)
    }
  }

  # 3d) invoke matchRanges()
  if (!is.null(seed)) set.seed(seed)
  mgr <- do.call(
    nullranges::matchRanges,
    c(
      list(
        focal = foreground_elements,
        pool  = background_universe,
        covar = covar_formula
      ),
      nullranges_args
    )
  )  #  [oai_citation:0‡bioconductor.org](https://www.bioconductor.org/packages//release/bioc/vignettes/nullranges/inst/doc/matchRanges.html)

  # 3e) extract the matched set as a GRanges
  matched_gr <- matched(mgr)  # same as pool(mgr)[indices(mgr, set="matched")]

  # 3f) if user wants a different ratio than 1:1, sample/truncate or pad
  n_fg <- length(foreground_elements)
  target_n <- n_ratio * n_fg
  current_n <- length(matched_gr)
  if (current_n > target_n) {
    matched_gr <- matched_gr[sample(current_n, target_n)]
  } else if (current_n < target_n && allow_replacement) {
    extra_idx <- sample(current_n, target_n - current_n, replace = TRUE)
    matched_gr <- c(matched_gr, matched_gr[extra_idx])
  }

  # 4) return the matched GRanges
  matched_gr
}
