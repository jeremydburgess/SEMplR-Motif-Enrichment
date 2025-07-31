#' Guess Likely ID Keytypes from User IDs
#'
#' @description
#' \code{helper_guessCols()} inspects a vector of user-supplied identifiers and
#' heuristically returns the most probable AnnotationDbi keytype(s) to use for mapping.
#'
#' @param ul Character vector of identifiers (e.g., gene or transcript IDs).
#'
#' @return Character vector of one or more keytype names. Possible returns:
#' \describe{
#'   \item{\dQuote{entrezid}}{if all entries are numeric; Entrez Gene IDs.}
#'   \item{\dQuote{ensembl, ensembltrans, ensemblprot}}{if all entries match ENSG/ENST/ENSP patterns.}
#'   \item{\dQuote{refseq}}{if all entries match N[MRP]_\d+ patterns (RefSeq IDs).}
#'   \item{\dQuote{symbol, alias}}{fallback for other ID formats.}
#' }
#'
#' @examples
#' helper_guessCols(c("1234", "5678"))                # returns "entrezid"
#' helper_guessCols(c("ENSG00000123456", "ENSG00000987654"))
#' # returns c("ensembl", "ensembltrans", "ensemblprot")
#'
#' @keywords internal
helper_guessCols <- function(ul) {
  if (all(grepl("^[0-9]+$", ul)))        return("entrezid")
  if (all(grepl("^ENS(G|T|P)\\d+", ul))) return(c("ensembl","ensembltrans","ensemblprot"))
  if (all(grepl("^(N[MRP]_)", ul)))      return("refseq")
  c("symbol","alias")
}

