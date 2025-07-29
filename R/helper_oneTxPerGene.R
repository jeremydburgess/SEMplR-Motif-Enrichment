#' Select one transcript per gene (widest, then 5′-most)
#'
#' From a `GRanges` with potentially multiple entries per gene, returns exactly
#' one range per unique `entrez` ID by:
#' 1. Selecting the maximum width.
#' 2. If there is a tie, choosing the range with the most upstream 5′ coordinate.
#'
#' @param gr A `GRanges` object with an `entrez` metadata column.
#' @return A `GRanges` object containing exactly one range per gene.
#' @importFrom dplyr mutate select group_by ungroup filter arrange slice_head
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export
helper_oneTxPerGene <- function(gr) {
 gr_df <- as.data.frame(gr)

 gr_df <- gr_df %>%
   dplyr::mutate(rangeStart = dplyr::if_else(strand == "+",start,end)) %>%
   dplyr::group_by(entrez) %>%
   dplyr::filter(width == max(width, na.rm = TRUE)) %>%
   dplyr::arrange(rangeStart) %>%
   dplyr::slice_head(n = 1) %>%
   dplyr::ungroup() %>%
   dplyr::select(-rangeStart)

 makeGRangesFromDataFrame(gr_df,keep.extra.columns = TRUE)

 }







