source("R/parseGENCODEgtf.R")
source("R/annotateGTF.R")
source("R/computeMatchStats.R")
source("R/mapUserList.R")
source("R/definePools.R")
source("R/definePromoterRegions.R")
source("R/defineBackgroundElements.R")
source("R/getPromoterSeqs.R")

library(readr)    # for read_tsv()
library(dplyr)    # for filter(), mutate(), distinct(), pull()
library(tidyr)    # for separate_rows()
library(stringr)  # for str_squish(), str_count(), str_match(), str_match_all(), word()




# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", version = "3.20")
# BiocManager::install("nullranges", version = "3.20")
library(BSgenome.Hsapiens.UCSC.hg38)
library(nullranges)
library(ks)


# Parse gencode gtf
gtf_path <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz"
gtf_path <- "/Users/jeremy_burgess/Downloads/gencode.v48.annotation.gtf.gz" # Locally stored version
gtf <- parseGENCODEgtf(feature_types = c("gene","transcript"))

# Add priority flags and alternate id types
# gtf_map <- annotateGTF(gtf)
gtf_map <- annotateGTF(gtf,metadata_urls = list( # Locally stored version
  entrez = "/Users/jeremy_burgess/Downloads/gencode.v48.metadata.EntrezGene.gz",
  hgnc   = "/Users/jeremy_burgess/Downloads/gencode.v48.metadata.HGNC.gz",
  refseq = "/Users/jeremy_burgess/Downloads/gencode.v48.metadata.RefSeq.gz"
))


# Based on provided input list, map to an id type and return minimal table to use to identify promoter coordinates for matching genes/transcripts

# Make dummy user list (500 random ids from gtf_map) - replace gene_id for other id types (doesn't account for NAs or forcing unique ids)
user_list <- gtf_map$gene_id[sample.int(nrow(gtf_map), 500, replace = FALSE)]

# Map user list to gtf_map and extract relevant columns from subset (foreground)
mappedRecords <- mapUserList(user_list,gtf_map)

# Filter foreground elements and background pool according to universe and gene filters
TSSforEnrichment <- definePools(mappedRecords,id_level = "auto",
                                gene_filter = "flag:Ensembl_canonical",
                                universe_filter = list(gene_type = c("lncRNA","protein_coding")),
                                universe_filter_logic = "or")



# Define promoter regions and convert to GRangesList (maintaining other columns as metadata)
PromoterGRanges <- definePromoterRegions(TSSforEnrichment)
    # PromoterGRanges$backgroundUniverse  # GRanges
    # PromoterGRanges$foregroundElements  # GRanges

# Select background elements from backgroundPool
## 1) Return the raw pool (including FG) ##
bg_pool <- defineBackgroundElements(
  background_universe = PromoterGRanges$backgroundUniverse,
  foreground_elements = PromoterGRanges$foregroundElements,
  bsgenome            = BSgenome.Hsapiens.UCSC.hg38,
  method              = "pool",
  exclude_fg          = FALSE
)


## 2) Random draw at 2× the FG size ##
bg_random <- defineBackgroundElements(
  background_universe = PromoterGRanges$backgroundUniverse,
  foreground_elements = PromoterGRanges$foregroundElements,
  bsgenome            = BSgenome.Hsapiens.UCSC.hg38,
  method              = "random",
  n_ratio             = 2,       # twice as many bg as fg
  exclude_fg          = TRUE,    # drop true promoters first
  seed                = 42
)


## 3) Matched draw on width + GC + chromosome ##
bg_matched <- defineBackgroundElements(
  background_universe = PromoterGRanges$backgroundUniverse,
  foreground_elements = PromoterGRanges$foregroundElements,
  bsgenome            = BSgenome.Hsapiens.UCSC.hg38,
  method              = "matched",
  n_ratio             = 1,
  covariates          = c("width","gc","seqnames"),  # match on length, GC, chr
  exclude_fg          = TRUE,
  seed                = 2025,
  nullranges_args     = list(
    method     = "rejection"
  )
)


# Add sequences to GRanges
PromoterGRangesSeqs <- getPromoterSeqs(matched_bg, BSgenome.Hsapiens.UCSC.hg38, flank_rc = TRUE)
