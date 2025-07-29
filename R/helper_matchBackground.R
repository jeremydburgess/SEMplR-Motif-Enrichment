#' Match Background Regions to Foreground by Key Covariates
#'
#' @description
#' Uses \pkg{nullranges} to select a 1:1 background set that matches the
#' foreground on specified covariates (e.g. binned width, GC content, chromosome).
#'
#' @note
#' When using `nrMethod = "rejection"`, continuous covariates like promoter width
#' often have very low variance (e.g. many regions exactly 350 bp), which can break
#' the kernel‐density–based sampler. To match on width, first **bin** it into
#' categories (e.g. quartiles) so it is treated as a factor. Otherwise, omit width
#' or use `nrMethod = "stratified"`.
#'
#' @param pool A \code{GRanges} of pruned background candidate regions.
#' @param focal A \code{GRanges} of foreground regions to match.
#' @param organism Character(1). Species name (e.g. "Homo sapiens"),
#'   used to resolve a BSgenome if `gc` is requested.
#' @param genomeBuild Character(1). UCSC‐style genome build (e.g. "hg38").
#' @param bsgenome A loaded BSgenome object (e.g. \code{BSgenome.Hsapiens.UCSC.hg38}).
#'   If NULL and `gc` is in `covariates`, the helper will attempt to resolve one.
#' @param covariates Character vector of columns to match on. Valid values:
#'   * "width" → injected as `promoterWidth`
#'   * "seqnames" or "chromosome" → injected as `chromosome`
#'   * "gc" → computed from `bsgenome`
#'   * any other existing `mcols(pool)` column provided by the user
#' @param nrMethod Character(1). Sampling algorithm passed to
#'   \code{\link[nullranges]{matchRanges}()}: \code{"rejection"} or \code{"stratified"}.
#' @param bgReplace Logical. Allow replacement when sampling background
#'   (defaults to FALSE). If \code{nrMethod = "nearest"}, replacement is forced TRUE.
#' @param seed Integer. Optional random seed for reproducibility (default NULL).
#'
#' @return A \code{MatchRanges} object, containing:
#'   \itemize{
#'     \item \code{matched(mgr)}: the matched background \code{GRanges}
#'     \item \code{focal(mgr)}: the foreground \code{GRanges}
#'     \item \code{pool(mgr)}: the pruned background universe
#'     \item \code{unmatched(mgr)}: any focal regions that failed to match
#'   }
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' mgr <- helper_matchBackground(
#'   pool       = prom_bg,
#'   focal      = prom_fg,
#'   organism   = "Homo sapiens",
#'   genomeBuild= "hg38",
#'   bsgenome   = BSgenome.Hsapiens.UCSC.hg38,
#'   covariates = c("width", "gc", "chromosome"),
#'   nrMethod   = "rejection",
#'   bgReplace  = FALSE,
#'   seed       = 42
#' )
#' # extract the GRangesList of fg+bg:
#' GRangesList(
#'   foregroundElements = nullranges::focal(mgr),
#'   backgroundElements = nullranges::matched(mgr)
#' )
#' }
#'
#' @importFrom nullranges matchRanges matched focal pool unmatched
#' @importFrom Biostrings getSeq letterFrequency
#' @importFrom GenomicRanges width seqnames
#' @importFrom S4Vectors mcols
#' @importFrom stats as.formula
#' @export
helper_matchBackground <- function(
    pool,
    focal,
    organism,
    genomeBuild,
    bsgenome,
    covariates, # character vector passed in by helper_defineBackgroundElements
    nrMethod, # character string passed in by helper_defineBackgroundElements
    bgReplace,
    seed
) {
  if (!is.null(seed)) set.seed(seed)

  # 1) Inject reserved covariates into mcols() and rename to avoid conflicts
  # width → promoterWidth
  if ("width" %in% covariates) {
    mcols(pool)$promoterWidth  <- width(pool)
    mcols(focal)$promoterWidth <- width(focal)
    covariates[covariates == "width"] <- "promoterWidth"
  }

  # If both seqnames & chromosome were requested, drop seqnames
  if (all(c("seqnames", "chromosome") %in% covariates)) {
    covariates <- setdiff(covariates, "seqnames")
  }

  # seqnames → chromosome
  if ("seqnames" %in% covariates) {
    mcols(pool)$chromosome  <- as.character(seqnames(pool))
    mcols(focal)$chromosome <- as.character(seqnames(focal))
    covariates[covariates == "seqnames"] <- "chromosome"
  }

  # chromosome → (reuse) chromosome  ← inject if not already
  if ("chromosome" %in% covariates) {
    if (!"chromosome" %in% names(mcols(pool))) {
      mcols(pool)$chromosome  <- as.character(seqnames(pool))
      mcols(focal)$chromosome <- as.character(seqnames(focal))
    }
  }

  # GC content
  if ("gc" %in% covariates) {
    if (is.null(bsgenome)) {
      bsgenome <- helper_resolveBSgenome(organism, genomeBuild)
    }
    seqs_bg <- getSeq(bsgenome, pool)
    seqs_fg <- getSeq(bsgenome, focal)
    mcols(pool)$gc  <- rowSums(letterFrequency(seqs_bg, "GC",   as.prob = TRUE))
    mcols(focal)$gc <- rowSums(letterFrequency(seqs_fg, "GC",   as.prob = TRUE))
  }

  # 2) Allow any extra user‐supplied mcols
  extra <- setdiff(covariates, c("width","gc","seqnames","chromosome"))
  if (length(extra) && !all(extra %in% names(mcols(pool)))) {
    stop("Unknown covariates: ", paste(setdiff(extra, names(mcols(pool))), collapse=", "))
  }

  # 3) Drop zero‐variance covariates
  keep <- vapply(covariates, function(col) {
    vals <- mcols(pool)[[col]]
    if (is.numeric(vals)) var(vals, na.rm=TRUE)>0 else length(unique(vals))>1
  }, logical(1))
  covariates <- covariates[keep]

  # 4) Build the formula and ensure argument compatibility
  covar_formula <- as.formula(paste0("~", paste(covariates, collapse=" + ")))

  if(nrMethod == "nearest" && bgReplace == FALSE){
    message("'nearest' method not available for 'bgReplace = FALSE'. Setting 'bgReplace' to 'TRUE'")
    bgReplace <- TRUE
  }

  # 5) Call matchRanges()
  mgr <- nullranges::matchRanges(
    focal   = focal,
    pool    = pool,
    covar   = covar_formula,
    method  = nrMethod,
    replace = bgReplace
  )

  return(mgr)
}
