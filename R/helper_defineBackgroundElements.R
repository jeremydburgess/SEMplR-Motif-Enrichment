#' Define and Sample Background Elements for Motif Enrichment
#'
#' @description
#' Filters and samples a set of background genomic ranges to match a given set
#' of foreground ranges based on user-specified covariates. Supports three methods:
#' * `pool`: returns the pruned universe without sampling
#' * `random`: draws a random sample of size `n_ratio * length(foreground_elements)`
#' * `matched`: uses `nullranges::matchRanges()` to perform 1:1 matching on covariates
#'
#' @param background_universe A `GRanges` of all candidate background regions (already pruned).
#' @param foreground_elements A `GRanges` of target regions for which backgrounds are needed.
#' @param bgMethod Character. One of `"pool"`, `"random"`, or `"matched"`.
#' @param n_ratio Numeric scalar. Number of background regions = `n_ratio * number of foreground elements`.
#' @param bgExcludeFgOverlaps Logical. If `TRUE`, removes any background ranges overlapping foreground.
#' @param bgExcludeFgGenes Logical. If `TRUE`, removes any background whose gene (via `mcols()$entrez`) is in the foreground.
#' @param covariates Character vector. Covariates to match on for `matched` method. Valid entries:
#'   * `"width"` &rarr; injected as `promoterWidth`
#'   * `"seqnames"` or `"chromosome"` &rarr; injected as `chromosome`
#'   * `"gc"` &rarr; computed from `bsGenome`
#'   * any existing `mcols(background_universe)` column
#' @param bgReplace Logical. Allow replacement when sampling (applies to `random` and `matched`).
#' @param seed Integer. Optional random seed for reproducibility.
#' @param nrMethod Character. Sampling algorithm for `matched`: one of `"rejection"`,`"nearest"`,`"stratified"`.
#' @param bsGenome A loaded `BSgenome` object (e.g., `BSgenome.Hsapiens.UCSC.hg38`).
#'   If `NULL` and `"gc"` is requested, the function will attempt to resolve one via `helper_resolveBSgenome()`.
#' @param organism Character(1). Species name (e.g., `"Homo sapiens"`), required if `bsGenome` must be resolved.
#' @param genomeBuild Character(1). UCSC build (e.g., `"hg38"`), required if resolving `bsGenome`.
#'
#' @return A named `list` with elements:
#' * `backgroundElements`: `GRanges` of sampled background regions
#' * `foregroundElements`: the input foreground `GRanges`
#' * `backgroundUniverse`: the pruned universe `GRanges`
#' * `matchObject`: a `MatchedGRanges` object if `bgMethod = "matched"`, else `NULL`
helper_defineBackgroundElements <- function(
    background_universe,
    foreground_elements,
    bgMethod    = c("pool", "random", "matched"),
    n_ratio             = 1,
    bgExcludeFgOverlaps     = TRUE,
    bgExcludeFgGenes        = FALSE,
    covariates      = c("width", "gc", "seqnames"),
    bgReplace          = FALSE,
    seed                = NULL,
    nrMethod = c("rejection","nearest","stratified"),
    bsGenome = NULL,
    organism,
    genomeBuild
) {
  bgMethod <- match.arg(bgMethod)
  nrMethod <- match.arg(nrMethod)

  ## 1) Pruning steps
  # Optional: remove all background ranges from pool if any overlap with foreground range
  if (bgExcludeFgOverlaps) {
    background_universe <- subsetByOverlaps(
      background_universe, foreground_elements, invert = TRUE
    )
  }

  # Optional: remove all background genes from pool if they appear in foreground
  if (bgExcludeFgGenes) {
    fg_genes <- unique(mcols(foreground_elements)$entrez)
    background_universe <- background_universe[
      ! mcols(background_universe)$entrez %in% fg_genes
    ]
  }

  selectedBg <- switch(
    bgMethod,
    pool = {background_universe},
    random = {helper_randomBackground(
      background_universe, foreground_elements,n_ratio,bgReplace,seed)},
    matched = {
      helper_matchBackground(background_universe,foreground_elements,organism,
                             genomeBuild,bsGenome,covariates,nrMethod,bgReplace,seed)}
  )

  # Standardize into bg_gr and capture full matchObject if present
  if (is(selectedBg, "MatchedGRanges")) {
    bg_gr       <- nullranges::matched(selectedBg)
    fg_gr       <- nullranges::focal(selectedBg)
    universe    <- nullranges::pool(selectedBg)
    matchObject <- selectedBg
  } else {
    bg_gr       <- selectedBg
    fg_gr       <- foreground_elements
    universe    <- background_universe
    matchObject <- NULL
  }

  # Return a consistent list object
  list(
    backgroundElements = bg_gr,
    foregroundElements = fg_gr,
    backgroundUniverse = universe,
    matchObject        = matchObject
  )

}






