source("R/parseGENCODEgtf.R")
source("R/annotateGTF.R")
source("R/computeMatchStats.R")
source("R/mapUserList.R")
source("R/assembleTSSforeground.R")
source("R/definePromoterRegions.R")

# Parse gencode gtf
gtf_path <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz"
gtf <- parseGENCODEgtf(feature_types = c("gene","transcript"))

# Add priority flags and alternate id types
gtf_map <- annotateGTF(gtf)

# Based on provided input list, map to an id type and return minimal table to use to identify promoter coordinates for matching genes/transcripts

# Make dummy user list (500 random ids from gtf_map) - replace gene_id for other id types (doesn't account for NAs or forcing unique ids)
user_list <- gtf_map$gene_id[sample.int(nrow(gtf_map), 500, replace = FALSE)]

# Map user list to gtf_map and extract relevant columns from subset (foreground)
mappedRecords <- mapUserList(user_list,gtf_map)

# Extract dataframe of mapped records
mapped <- mappedRecords$mapped_df

# Filter to enrichment foreground TSS table
ForegroundTSS <- assembleTSSforeground(mappedRecords,id_level = "auto",gene_filter = "flag:Ensembl_canonical")

# Define promoter regions and convert to GRanges (maintaining other columns as metadata)
ForegroundPromoters <- definePromoterRegions(ForegroundTSS)
