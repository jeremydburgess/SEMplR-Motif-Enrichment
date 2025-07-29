#' Resolve and load a BSgenome for a given UCSC build
#'
#' Looks first among installed BSgenome data packages, then falls back to
#' Bioconductor’s available genomes list (and installs if needed).
#'
#' @param genomeBuild Character(1). UCSC‐style genome build (e.g. `"hg38"`).
#' @param bsMasking Logical(1). Whether to require the “masked” version
#'   (i.e. masked repeats/ambiguous regions) (default `FALSE`).
#' @param suffix Optional character. A package name suffix (e.g. `"dbSNP151.major"`),
#'   to disambiguate multiple hg38 variants.
#' @return A loaded \pkg{BSgenome} object.
#' @examples
#' \dontrun{
#'   # load your installed hg38 BSgenome or install it if needed
#'   bs <- helper_resolveBSgenome("hg38")
#'   # get the masked version if you prefer
#'   bs_masked <- helper_resolveBSgenome("hg38", bsMasking = TRUE)
#'   # or explicitly ask for dbSNP151.major
#'   bs_major <- helper_resolveBSgenome("hg38", suffix = "dbSNP151.major")
#' }
#' @importFrom BSgenome installed.genomes available.genomes
#' @importFrom BiocManager install
#' @export
helper_resolveBSgenome <- function(organism, genomeBuild, bsMasking=FALSE, suffix=NULL) {
  # 1) Look in installed.genomes()
  inst <- BSgenome::installed.genomes(splitNameParts=TRUE)
  cand <- inst[inst$organism == helper_bsOrganismCode(organism) &
                 inst$genome == genomeBuild &
                 inst$masked == bsMasking, ]
  if (nrow(cand) > 0) {
    pkg <- cand$pkgname[1]
    bsgenome_obj <- get(pkg, envir = asNamespace(pkg))
    return(bsgenome_obj)
  }

  # 2) Otherwise look in available.genomes()
  avail <- BSgenome::available.genomes(splitNameParts=TRUE)
  cand <- avail[avail$organism == helper_bsOrganismCode(organism) &
                  avail$genome == genomeBuild &
                  avail$masked == bsMasking, ]
  if (!is.null(suffix)) {
    cand <- subset(cand, grepl(suffix, pkgname))
  }
  if (nrow(cand) == 0) {
    stop("No BSgenome found for ", organism, " + ", genomeBuild)
  }
  pkg <- cand$pkgname[1]

  # 3) Install & load
  BiocManager::install(pkg)
  # After installation:
  requireNamespace(pkg, quietly = TRUE)
  bsgenome_obj <- get(pkg, envir = asNamespace(pkg))
  return(bsgenome_obj)
  }


