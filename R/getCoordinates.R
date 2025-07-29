#' Retrieve Promoter Regions and Sample Background Elements
#'
#' @description
#' Given a mapped set of foreground IDs (with Entrez and mappedID),
#' this function:
#' 1. Fetches genomic coordinates for each feature (gene or transcript), according
#'    to the chosen \code{TSS.method}.
#' 2. Extracts promoter windows around those coordinates.
#' 3. Optionally reduces overlaps and selects one promoter per gene.
#' 4. Samples background promoter elements matching the foreground, via
#'    \code{\link{helper_defineBackgroundElements}()}.
#'
#' @param mapped A list as returned by \code{\link{mapForegroundIDs}()}, containing
#'   at least \code{fg_ids}, \code{bg_ids}, \code{so_obj}, \code{orgdb},
#'   \code{genomeBuild}, and \code{organism}.
#' @param transcript Logical; \code{TRUE} for transcript‐level coordinates,
#'   \code{FALSE} for gene‐level.
#' @param TSS.method Character(1). Method to pick TSS when \code{transcript=FALSE}:
#'   \code{"UCSCgene"}, \code{"Ensembl_canonical"}, \code{"commonTSS"},
#'   \code{"uniqueTSS"}, \code{"fivePrimeTSS"}, or \code{"allTSS"}.
#' @param bgMethod Character; one of \code{"pool"}, \code{"random"}, or
#'   \code{"matched"}. Passed to \code{\link{helper_defineBackgroundElements}()}.
#' @param n_ratio Numeric; for \code{bgMethod="random"}, number of backgrounds =
#'   \code{n_ratio * length(foreground)}.
#' @param bgExcludeFgOverlaps Logical; if \code{TRUE}, drop any background whose
#'   promoter overlaps a foreground.
#' @param bgExcludeFgGenes Logical; if \code{TRUE}, drop any background gene
#'   whose EntrezID is in the foreground.
#' @param covariates Character vector of covariate names for matching:
#'   \code{"width"}, \code{"gc"}, \code{"seqnames"} (or \code{"chromosome"}), or
#'   any user‐provided \code{mcols()} column.
#' @param bgReplace Logical; allow replacement when sampling backgrounds.
#' @param seed Integer; optional random seed.
#' @param nrMethod Character; \code{"rejection"}, \code{"nearest"}, or
#'   \code{"stratified"}. Passed to \code{\link{helper_defineBackgroundElements}()}.
#' @param bsGenome A \pkg{BSgenome} object. If \code{NULL} and \code{"gc"} is in
#'   \code{covariates}, will be auto‐resolved via \code{\link{helper_resolveBSgenome}()}.
#' @param promoterWindow Numeric named vector of lengths:
#'   \code{c(upstream, downstream)} (default \code{c(300,50)}).
#' @param standardChroms Logical; restrict to standard chromosomes.
#' @param reduceOverlaps Logical; merge any overlapping promoter windows.
#' @param overlapMinGap Numeric; minimum gap when reducing overlaps.
#' @param onePromoterPerGene Logical; if \code{TRUE}, choose one promoter per gene.
#' @param ensDb Optional \code{EnsDb} object; required for
#'   \code{TSS.method="Ensembl_canonical"}.
#'
#' @return A named \code{list} from the final \code{helper_defineBackgroundElements()}:
#'   \item{backgroundElements}{\code{GRanges} of sampled background promoters}
#'   \item{foregroundElements}{\code{GRanges} of foreground promoters}
#'   \item{backgroundUniverse}{\code{GRanges} of the pruned universe}
#'   \item{matchObject}{\code{MatchedGRanges} if \code{bgMethod="matched"}, else \code{NULL}}
#'
#' @examples
#' \dontrun{
#' mapping <- buildMappingObject("Homo sapiens")
#' mapped  <- mapForegroundIDs(my_ids, mapping)
#' filtered <- poolFilter(mapped, geneType="protein-coding")
#' coords <- getCoordinates(
#'   filtered,
#'   transcript     = FALSE,
#'   TSS.method     = "UCSCgene",
#'   bgMethod       = "matched",
#'   n_ratio        = 1,
#'   covariates     = c("width","gc","chromosome"),
#'   nrMethod       = "stratified",
#'   bsGenome       = BSgenome.Hsapiens.UCSC.hg38,
#'   promoterWindow = c(upstream=300, downstream=50)
#' )
#' }
#'
#' @importFrom dplyr tbl select filter collect inner_join transmute group_by slice_head ungroup
#' @importFrom GenomicRanges makeGRangesFromDataFrame promoters
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom AnnotationDbi select
#' @importFrom S4Vectors mcols
#' @export
getCoordinates <- function(mapped,
                           transcript = FALSE,
                           TSS.method = c("UCSCgene","Ensembl_canonical",
                                          "commonTSS","uniqueTSS","fivePrimeTSS","allTSS"),
                           bgMethod             = c("pool", "random", "matched"),
                           n_ratio              = 1,
                           bgExcludeFgOverlaps  = TRUE,
                           bgExcludeFgGenes     = FALSE,
                           covariates           = c("width", "gc", "chromosome"),
                           bgReplace            = FALSE,
                           seed                 = NULL,
                           nrMethod             = c("rejection","nearest","stratified"),
                           bsGenome             = NULL,
                           promoterWindow       = c(upstream=300, downstream=50),
                           standardChroms       = TRUE,
                           reduceOverlaps       = TRUE,
                           overlapMinGap        = 0,
                           onePromoterPerGene   = TRUE,
                           ensDb                = NULL) {

  TSS.method <- match.arg(TSS.method)
  bgMethod <- match.arg(bgMethod)
  nrMethod <- match.arg(nrMethod)

  # unpack
  fg_df       <- mapped$fg_ids      # data.frame(entrez, mappedID)
  bg_df       <- mapped$bg_ids
  so_obj      <- mapped$so_obj
  transcript  <- mapped$transcript   # logical flag
  genomeBuild <- mapped$genomeBuild
  organism    <- mapped$organism

  # first make sure we have an EnsDb if the user asked for Ensembl_canonical
  if (TSS.method == "Ensembl_canonical") {
    if (is.null(ensDb) && !is.null(mapped$ensdb)) {
      ensDb <- mapped$ensdb
    }
    if (is.null(ensDb)) {
      # attempt to auto-load one via your helper
      ensDb <- tryCatch(helper_loadEnsDbForOrganism(mapped$organism, mapped$genomeBuild),
        error = function(err) {
          stop("TSS.method = 'Ensembl_canonical' requires an EnsDb object.\n",
            "No suitable EnsDb found for ", mapped$organism," + ", mapped$genomeBuild, ".") }
        )} }

  if (transcript) {
    #
    # === TRANSCRIPT‐LEVEL MODE ===
    # We know for it to have got this far, transcript = TRUE means mappedID is ensembltrans
    # So we can directly reference transcript tables with no extra mapping required

    # Get relevant "ranges_tx" table with version stripped ensembltrans (tx_name)
    ranges_tx <- tbl(so_obj, "ranges_tx") %>%
            dplyr::select(tx_name,tx_chrom,tx_start,tx_end,tx_strand) %>%
      dplyr::collect () %>%
      dplyr::transmute(ensembltrans = sub("\\.\\d+$", "", tx_name),  # strip off “.1” suffix *and* rename
        tx_chrom,tx_start,tx_end,tx_strand)

    # Join bg_df with ranges_tx on ensembltrans
    bg_ranges <- dplyr::inner_join(bg_df,ranges_tx, # inner rather than left as no point keeping rows with no coords
                               by = dplyr::join_by(mappedID == ensembltrans)) %>% unique()

    # build GRanges for bg pool
    gr_bg <- GenomicRanges::makeGRangesFromDataFrame(
      bg_ranges,
      seqnames.field     = "tx_chrom",
      start.field        = "tx_start",
      end.field          = "tx_end",
      strand.field       = "tx_strand",
      keep.extra.columns = TRUE
    )

    # Restrict to standard chromosomes? Optional but default yes
    if(standardChroms){
      gr_bg <- gr_bg %>%
        GenomeInfoDb::keepStandardChromosomes(species = stringr::str_to_sentence(
                        stringr::str_to_lower(organism)),pruning.mode = "coarse")
    }
    # Generate foreground elements granges by subsetting bg_gr by mappedID
    gr_fg <- gr_bg[gr_bg$mappedID %in% fg_df$mappedID]

    } else {

  #
  # === GENE‐LEVEL MODE ===
  #

  # Deal with UCSCgene first since it requires joining to a different table
      if (TSS.method == "UCSCgene") {
    # get relevant gene_ranges table
    ranges_gene <- tbl(so_obj, "ranges_gene") %>%
      dplyr::select(entrez,gene_chrom,gene_start,gene_end,gene_strand) %>%
      dplyr::collect()

    # if standardChrom = TRUE (default behavior), restrict to standard chroms
    if(standardChroms) {ranges_gene <- dplyr::filter(ranges_gene, gene_chrom %in%
                                                     subset(getChromInfoFromUCSC(genomeBuild),
                                                            assembled == TRUE)$chrom)}

    # Join bg_df with ranges_gene on entrez
    bg_ranges_gene <- dplyr::inner_join(bg_df,ranges_gene, # inner rather than left as no point keeping rows with no coords
                                   by = "entrez", relationship = "many-to-many") %>% unique()

    # Build GRanges for bg
    gr_bg <- GenomicRanges::makeGRangesFromDataFrame(
      bg_ranges_gene,
      seqnames.field     = "gene_chrom",
      start.field        = "gene_start",
      end.field          = "gene_end",
      strand.field       = "gene_strand",
      keep.extra.columns = TRUE
    )

    # Generate foreground elements granges by subsetting bg_gr by mappedID
    gr_fg <- gr_bg[gr_bg$entrez %in% fg_df$entrez]

    } else {

  # for all other methods (they require transcript‐level TSS computations)
  # We'll map entrez from bg_df to entrez in ranges_tx src table
  # First get relevant columns from ranges_tx src table
  ranges_tx <- tbl(so_obj, "ranges_tx") %>%
    dplyr::select(entrez,tx_name,tx_chrom,tx_start,tx_end,tx_strand) %>%
    dplyr::collect () %>% unique()

  # Join bg_df with ranges_tx on ensembltrans
  bg_ranges <- dplyr::inner_join(bg_df,ranges_tx, # inner rather than left as no point keeping rows with no coords
                                 by = "entrez",relationship = "many-to-many") %>% unique()

# Compute TSS based on strand
  bg_ranges$tss <- ifelse(bg_ranges$tx_strand == "+", bg_ranges$tx_start, bg_ranges$tx_end)

  # if standardChrom = TRUE (default behavior), restrict to standard chroms
  if(standardChroms) {bg_ranges <- dplyr::filter(bg_ranges, tx_chrom %in%
                                                   subset(getChromInfoFromUCSC(genomeBuild),
                                                          assembled == TRUE)$chrom)}

  # Now perform specific TSS-based behavior based on TSS.method
  # NB no tie-breaking here, and multiple tss per gene are allowed
  tss_df <- switch(
    TSS.method,
    Ensembl_canonical = helper_filterByEnsemblCanonical(bg_ranges,ensDb),
    commonTSS   = helper_selectCommonTSS(bg_ranges),
    uniqueTSS   = bg_ranges %>% dplyr::group_by(entrez,tss) %>% dplyr::slice_head(n=1) %>% ungroup(),
    fivePrimeTSS = helper_selectFivePrimeTSS(bg_ranges), # if tied, just picks a single record
    allTSS      = bg_ranges
  )

  # 5) build bg GRanges
gr_bg <- GenomicRanges::makeGRangesFromDataFrame(
    tss_df,
    seqnames.field     = "tx_chrom",
    start.field        = "tx_start",
    end.field          = "tx_end",
    strand.field       = "tx_strand",
    keep.extra.columns = TRUE
  )

  # filter to fg from fg_df
  gr_fg <- gr_bg[gr_bg$entrez %in% fg_df$entrez]
}
      }

  ####── promoter extraction ──────────────────────────────────────────####
  # gr_bg, gr_fg are GRanges with metadata column “entrez” (gene ID)

  upstream = promoterWindow[["upstream"]]
  downstream = promoterWindow[["downstream"]]

  prom_bg <- GenomicRanges::promoters(
    gr_bg,
    upstream   = upstream,
    downstream = downstream,
    use.names  = TRUE
  )
  prom_fg <- GenomicRanges::promoters(
    gr_fg,
    upstream   = upstream,
    downstream = downstream,
    use.names  = TRUE
  )


# (optional) merge any overlapping promoter windows within each set
  if (reduceOverlaps) {
    prom_bg <- helper_reduceOverlapsWithinGenes(prom_bg,overlapMinGap)
    prom_fg <- helper_reduceOverlapsWithinGenes(prom_fg,overlapMinGap)
  }

  # (optional) enforce one promoter per gene: choose the widest window, break ties with 5'-most
  if (onePromoterPerGene) {
    prom_bg <- helper_oneTxPerGene(prom_bg)
    prom_fg <- helper_oneTxPerGene(prom_fg)
  }


  ####── background element selection ──────────────────────────────────────────####

  out <- helper_defineBackgroundElements(background_universe = prom_bg,
                                  foreground_elements = prom_fg,
                                  bgMethod            = bgMethod,
                                  n_ratio             = n_ratio,
                                  bgExcludeFgOverlaps = bgExcludeFgOverlaps,
                                  bgExcludeFgGenes    = bgExcludeFgGenes,
                                  covariates          = covariates,
                                  bgReplace           = bgReplace,
                                  seed                = seed,
                                  nrMethod            = nrMethod,
                                  bsGenome            = bsGenome)

return(out)

  }



