#' Filter a transcript‐coordinate table to Ensembl canonical transcripts
#'
#' @description
#' \code{helper_filterByEnsemblCanonical()} takes a data frame of transcript
#' coordinates (with a \code{tx_name} column), queries an \code{EnsDb} for
#' the set of canonical transcript IDs, and returns only those rows whose
#' (version‐stripped) \code{tx_name} matches a canonical transcript.
#'
#' @param df     A data.frame or tibble containing at least a column
#'               \code{tx_name}, which holds versioned Ensembl transcript IDs
#'               (e.g. “ENST0000…​.1”).
#' @param ensdb  An \code{EnsDb} object (from AnnotationHub or elsewhere).
#'
#' @return A data.frame containing only the rows of \code{df} whose transcript
#'         (after stripping the “.\d+” suffix) is marked canonical in \code{ensdb}.
#'
#' @details
#' Internally, this helper:
#' 1. Uses \code{AnnotationDbi::select()} on \code{ensdb} to retrieve all
#'    \code{TXID} ↔ \code{TXISCANONICAL} mappings.
#' 2. Filters to those with \code{TXISCANONICAL == 1}.
#' 3. Strips the trailing “.\d+” from \code{tx_name} in \code{df} to get base
#'    transcript IDs.
#' 4. Keeps only those rows whose base IDs appear in the canonical set.
#'
#' @importFrom AnnotationDbi select keys
#' @importFrom dplyr select mutate filter distinct pull
#' @keywords internal
helper_filterByEnsemblCanonical <- function(df,ensdb) {

# 1) get canonical transcript IDs from EnsDb
canonicalTranscripts <- AnnotationDbi::select(
  ensdb,
  keys = keys(ensdb,keytype = "TXID"),
  keytype = "TXID",
  columns = c("TXID","TXISCANONICAL")
)  %>% dplyr::filter(TXISCANONICAL == 1) %>%
  dplyr::select(TXID) %>%
  dplyr::pull()

# 2) strip version suffix and filter
df <- dplyr::mutate(df,tx_id = sub("\\.\\d+$", "", tx_name)) %>%
  dplyr::filter(tx_id %in% canonicalTranscripts) %>% unique()

return(df)
}
