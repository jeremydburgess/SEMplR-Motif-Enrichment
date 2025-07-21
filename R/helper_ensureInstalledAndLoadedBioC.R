#' Ensure a Bioconductor package is installed and loaded
#'
#' @description
#' Checks whether the specified Bioconductor package is installed; if not,
#' installs it via \code{BiocManager::install()}, then loads it into the
#' session with \code{library()}.
#'
#' @param pkg Character scalar giving the name of a Bioconductor package.
#'
#' @return
#' Invisibly returns \code{NULL}.  Called for its side effects of
#' installing and loading the package.
#'
#' @examples
#' \dontrun{
#' # Will install and load TxDb.Hsapiens.UCSC.hg38.knownGene if not already available
#' helper_ensureInstalledAndLoadedBioC("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' }
#'
#' @keywords internal
#' @importFrom BiocManager install
#' @export
helper_ensureInstalledAndLoadedBioC <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing Bioconductor package ‘", pkg, "’ …")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  # now load it
  library(pkg, character.only = TRUE)
}
