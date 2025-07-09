library(nullranges)    # for matchRanges()
library(Biostrings)    # for getSeq(), letterFrequency()
library(GenomicRanges)
library(ks) # for rejection method in matchRanges

#' Define Background Elements for Motif Enrichment
#'
#' @param background_universe A GRanges of all candidate regions.
#' @param foreground_elements A GRanges of the target (foreground) regions.
#' @param bsgenome A loaded BSgenome object, e.g. \code{BSgenome.Hsapiens.UCSC.hg38}.
#' @param method One of \code{"pool"}, \code{"random"}, or \code{"matched"}.
#' @param ...
#' @return GRangesList with \code{backgroundElements} and \code{foregroundElements}.
#' @export
defineBackgroundElements <- function(
    background_universe,    # GRanges of all candidate regions
    foreground_elements,    # GRanges of your true promoters
    bsgenome,                # BSgenome object, e.g. BSgenome.Hsapiens.UCSC.hg38"
    method            = c("pool","random","matched"),
    n_ratio           = 1,                      # #bg = n_ratio * #fg
    covariates        = c("width","gc","seqnames"),        # covariates to match (or existing mcols)
    exclude_fg        = TRUE,                   # drop fg from universe before sampling
    allow_replacement = FALSE,                  # for sampling
    seed              = NULL,                   # for reproducibility
    nullranges_args   = list(                   # passed to matchRanges()
      method     = "rejection",
      replace    = FALSE
    )
) {
  method <- match.arg(method)

  ## 1) Prune universe if requested
  if (exclude_fg) {
    background_universe <- subsetByOverlaps(
      background_universe, foreground_elements, invert = TRUE
    )
  }

  ## 2) Dispatch on method
  if (method == "pool") {
    bg_gr <- background_universe

  } else if (method == "random") {
    if (!is.null(seed)) set.seed(seed)
    n_bg <- n_ratio * length(foreground_elements)
    idx <- sample(
      seq_along(background_universe),
      size    = n_bg,
      replace = allow_replacement
    )
    bg_gr <- background_universe[idx]

  } else if (method == "matched") {
    # enforce 1:1 for true matching
    if (n_ratio != 1) {
      warning(
        "matchRanges() only returns a 1:1 background:focal set; ",
        "ignoring n_ratio = ", n_ratio, " and using n_ratio = 1 instead."
      )
      n_ratio <- 1
    }

    if (!requireNamespace("ks", quietly=TRUE)) {
      stop("Please install the 'ks' package to use matched sampling: install.packages('ks')")
    }

    # before injection, remap 'seqnames' to 'chr' (seqnames is reserved), width to featureWidth for same reason
    covariates[covariates == "seqnames"] <- "chr"
    covariates[covariates == "width"] <- "featureWidth"

    ## --- compute or inject each requested covariate ---
    for (cov in covariates) {
      if (cov == "featureWidth") {
        if (!"featureWidth" %in% names(mcols(background_universe))) {
          mcols(background_universe)$featureWidth        <- width(background_universe)
          mcols(foreground_elements)$featureWidth        <- width(foreground_elements)
        }

      } else if (cov == "gc") {
        if (!"gc" %in% names(mcols(background_universe))) {
          seqs_bg <- getSeq(bsgenome,         background_universe)
          seqs_fg <- getSeq(bsgenome,         foreground_elements)
          mcols(background_universe)$gc        <- rowSums(letterFrequency(seqs_bg, "GC",   as.prob=TRUE))
          mcols(foreground_elements)$gc        <- rowSums(letterFrequency(seqs_fg, "GC",   as.prob=TRUE))
        }

      } else if (cov == "chr") {
        if (!"chr" %in% names(mcols(background_universe))) {
          mcols(background_universe)$chr  <- as.character(seqnames(background_universe))
          mcols(foreground_elements)$chr  <- as.character(seqnames(foreground_elements))
        }

      } else if (cov %in% names(mcols(background_universe))) {
        # already present, do nothing

      } else {
        stop("Unknown covariate '", cov,
             "'. Must be 'featureWidth', 'gc', 'seqnames', 'chr', or an existing metadata column.")
      }
    }

    # coerce any character covariates into factors
    for (cov in covariates) {
      if (is.character(mcols(background_universe)[[cov]])) {
        mcols(background_universe)[[cov]] <- factor(mcols(background_universe)[[cov]])
        mcols(foreground_elements)[[cov]] <- factor(mcols(foreground_elements)[[cov]])
      }
    }

    # compute per‐covariate variance (or unique count)
    drop <- vapply(covariates, function(col) {
      vals <- mcols(background_universe)[[col]]
      if (is.numeric(vals)) var(vals)==0 else length(unique(vals))==1
    }, logical(1))
    if (any(drop)) {
      message("Dropping zero‐variance covariates: ", paste(covariates[drop], collapse=", "))
      covariates <- covariates[!drop]
    }

    # build the matching formula
    covar_formula <- as.formula(paste0("~", paste(covariates, collapse=" + ")))

    # filter out unsupported nullranges args
    valid_args <- names(formals(nullranges::matchRanges))
    nr_args    <- nullranges_args[names(nullranges_args) %in% valid_args]

    # perform the matchRanges call
    if (!is.null(seed)) set.seed(seed)
    mgr <- do.call(
      nullranges::matchRanges,
      c(
        list(
          focal = foreground_elements,
          pool  = background_universe,
          covar = covar_formula
        ),
        nr_args
      )
    )

    # extract matched set
    matched_gr <- matched(mgr)
    bg_gr <- matched_gr

  } else {
    stop("Unknown method: ", method)
  }

  ## 3) Return a GRangesList bundling BG + FG
  GRangesList(
    backgroundElements = bg_gr,
    foregroundElements = foreground_elements
  )
}
