

#### Test the pipeline

### Parse GENCODE file for mapping and coordinates
source("parseGENCODEgtf.R")

gtf_path <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz"
gtf <- parseGENCODEgtf(gtf_path,feature_types = c("gene","transcript"))


### Add tag flags (GENCODE_Primary, MANE_Select, Ensembl_canonical, MANE_Plus_Clinical)
gtf_flags <- gtf %>%
  mutate(GENCODE_Primary = ifelse(!is.na(tag) & str_detect(tag, "GENCODE_Primary"),1,0),
         Ensembl_canonical = ifelse(!is.na(tag) & str_detect(tag, "Ensembl_canonical"),1,0),
         MANE_Select = ifelse(!is.na(tag) & str_detect(tag, "MANE_Select"),1,0),
         MANE_Plus_Clinical = ifelse(!is.na(tag) & str_detect(tag, "MANE_Plus_Clinical"),1,0)
  )


### Add additional ID types
# Generate base ensembl gene and transcript IDs
gtf_extID <- gtf_flags %>%
  mutate(gene_id_base = sub("\\.\\d+$", "", gene_id),
    transcript_id_base = sub("\\.\\d+$", "", transcript_id),
    unique_id = row_number()
    )

# Import Gencode metadata for refseq, entrez, and hgnc
entrez_meta_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.metadata.EntrezGene.gz"
HGNC_meta_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.metadata.HGNC.gz"
RefSeq_meta_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.metadata.RefSeq.gz"

entrez_meta <- read_tsv(entrez_meta_url, col_names = c("transcript_id", "entrez_id"), show_col_types = FALSE)
HGNC_meta <- read_tsv(HGNC_meta_url, col_names = c("transcript_id", "HGNC_symbol", "HGNC_id"), show_col_types = FALSE)
RefSeq_meta <- read_tsv(RefSeq_meta_url, col_names = c("transcript_id", "RefSeq_TxId", "RefSeq_PxId"), show_col_types = FALSE)

# Join refseq, entrez, and hgnc
gtf_map <- gtf_extID %>%
  left_join(entrez_meta, by = "transcript_id") %>%
  left_join(HGNC_meta, by = "transcript_id") %>%
  left_join(RefSeq_meta, by = "transcript_id") %>%
  distinct()

# Identify supported ID columns
# Support: gene_name & transcript_name (What is the gene name in the GTF/GFF3?
                  #  Gene names are usually HGNC or MGI-approved gene
                  #  symbols mapped to the GENCODE genes by the Ensembl xref pipeline.
                  #  Sometimes, when there is no official gene symbol,
                  #  the Havana clone-based name is used.)
#           transcript_id (ensembl, versioned and _base)
#           gene_id (ensembl, versioned and _base)
#           entrez_id
#           HGNC (_symbol and _id)
#           RefSeq (_TxID and _PxID)

idCols <- c("unique_id", #Arbitrary row number to anchor selected id types to
            "gene_name",
            "gene_id",
            "gene_id_base",
            "entrez_id",
            "HGNC_symbol",
            "HGNC_id",
            "transcript_name",
            "transcript_id",
            "transcript_id_base",
            "RefSeq_TxID",
            "RefSeq_PxID")


# For each ID column, count number of InputIds and calculate proportion that exactly match.
compute_match_stats <- function(user_list,
                                map_df,
                                idCols = NULL {
  if (is.null(idCols)) {
    idCols <- c(
      "unique_id", "gene_name", "gene_id", "gene_id_base",
      "entrez_id", "HGNC_symbol", "HGNC_id", "transcript_name",
      "transcript_id", "transcript_id_base", "RefSeq_TxID", "RefSeq_PxID"
    )
  }
  # Only keep ID columns that actually exist
  cols <- intersect(idCols, names(map_df))
  total_ids <- length(user_list)

  stats <- lapply(cols, function(col) {
    # Unique non-NA values in this column
    vals <- unique(map_df[[col]])
    # Count how many user IDs matched
    n_matched <- sum(user_list %in% vals, na.rm = TRUE)
    # Compute percent
    pct <- n_matched / total_ids * 100

    data.frame(
      id_type      = col,
      matched      = n_matched,
      total        = total_ids,
      pct_matched  = pct,
      match_str    = paste0(n_matched, " of ", total_ids),
      stringsAsFactors = FALSE
    )
  })

  # Combine into one data.frame
  do.call(rbind, stats)
}

# Now identify which id type had the best mapping, or error if none are sufficient
map_and_report <- function(user_list, map_df, idCols, threshold = 0.90) {
  # 1) Compute match stats
  stats <- compute_match_stats(user_list, map_df, idCols)
  # stats <- compute_match_stats(user_list, gtf_map, idCols)

  # 2) Pick the best‐matching ID type
  best <- stats[which.max(stats$pct_matched), ]

  # 3) Check threshold
  if (best$pct_matched / 100 >= threshold) {
    msg <- sprintf(
      "Identified input ids as %s (%d of %d; %.1f%%).",
      best$id_type, best$matched, best$total, best$pct_matched
    )
    message(msg)

  # 4) Classify the ID type
  gene_keys       <- c("gene_id", "gene_name", "entrez_id", "HGNC_symbol", "HGNC_id")
  transcript_keys <- c("transcript_id", "transcript_name", "RefSeq_TxID")

  id_level <- if (best$id_type %in% gene_keys) {
    "gene"
  } else if (best$id_type %in% transcript_keys) {
    "transcript"
  } else {
    # fallback: assume transcript-level for anything else
    "transcript"
  }

  # 5) Subset to only that feature_type
  subset_df <- map_df %>%
    filter(feature_type == id_level,
           .data[[best$id_type]] %in% user_list)

  # 6) Build your minimal table
  minimal <- subset_df %>%
    select(
      unique_id, seqnames,
      start, end, strand
    ) %>%
    unique()

  return(minimal)
  } else {
    message(
      sprintf(
        "Unable to map IDs ≥ %.0f%% in any supported field; please supply IDs in one of: %s.",
        threshold*100,
        paste(stats$id_type, collapse = ", ")
      )
    )
    return(NULL)
  }
}




# Example usage:
user_list <- c("BRCA1","TP53","ENST00000367770","FOO")
# user_list <- c(user_list,gtf_map$gene_name[1:50])
user_list <- c(user_list,gtf_map$transcript_id[1:50])

matched <- map_and_report(user_list, gtf_map, idCols)



user_list











# Report back best choice of mapping with percentage, or if none above 90% say it can't be matched.
