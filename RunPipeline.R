source("R/parseGENCODEgtf.R")
source("R/annotateGTF.R")
source("R/computeMatchStats.R")
source("R/mapUserList.R")
# source("R/assembleTSSforeground.R")
source("R/definePools.R")
source("R/definePromoterRegions.R")
source("R/getPromoterSeqs.R")

# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", version = "3.20")
library(BSgenome.Hsapiens.UCSC.hg38)


# Parse gencode gtf
# gtf_path <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz"
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

# # # Extract dataframe of mapped records
# mapped <- mappedRecords$mapped_df
# universe <- mappedRecords$input_df

# # Filter to enrichment foreground TSS table
# ForegroundTSS <- assembleTSSforeground(mappedRecords,id_level = "auto",gene_filter = "flag:Ensembl_canonical")

# Filter foreground elements and background pool according to universe and gene filters
TSSforEnrichment <- definePools(mappedRecords,id_level = "auto",
                                gene_filter = "flag:Ensembl_canonical",
                                universe_filter = list(gene_type = c("lncRNA","protein_coding")),
                                universe_filter_logic = "or")



# Define promoter regions and convert to GRangesList (maintaining other columns as metadata)
PromoterGRanges <- definePromoterRegions(TSSforEnrichment)
    # PromoterGRanges$backgroundUniverse  # GRanges
    # PromoterGRanges$foregroundElements  # GRanges

# Add sequences to GRanges
PromoterGRangesSeqs <- getPromoterSeqs(PromoterGRanges, BSgenome.Hsapiens.UCSC.hg38, flank_rc = TRUE)








# Testing:

library(GenomicRanges)
library(tibble)

# say fg_prom is your GRanges and prom_seqs its DNAStringSet
df_seqs <- as_tibble(ForegroundPromoters) %>%       # turns GRanges → tibble with seqnames, start, end, strand, plus mcols
  mutate(sequence = as.character(prom_seqs))

df_seqs
