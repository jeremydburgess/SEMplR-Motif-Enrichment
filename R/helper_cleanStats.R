#' Clean up mapping statistics data.frame
#'
#' @description
#' Takes a data.frame with columns \code{matched} (integer counts per keytype)
#' and computes a new column \code{pct_matched} as the percentage of the user’s
#' input IDs that were mapped.  Any \code{NA} values in \code{matched} or
#' \code{pct_matched} are set to zero.
#'
#' @param df A data.frame containing at least:
#'   \describe{
#'     \item{\code{matched}}{Integer count of mapped IDs per keytype.}
#'   }
#'
#' @return
#' The same data.frame, with two modifications:
#'   \itemize{
#'     \item A new column \code{pct_matched} = \code{matched} / total_inputs * 100.
#'     \item Any \code{NA} in \code{matched} or \code{pct_matched} replaced by 0.
#'   }
#'
#' @details
#' The function uses the length of the original \code{foreground_ids} vector
#' (from the parent scope) to compute percentages.  Ensure that
#' \code{foreground_ids} is in scope when calling.
#'
#' @keywords internal
#' @export
helper_cleanStats <- function(df) {
  df$pct_matched <- df$matched / length(foreground_ids) * 100
  df$pct_matched[is.na(df$pct_matched)] <- 0
  df$matched[is.na(df$matched)] <- 0
  df
}
