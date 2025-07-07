# ──────────────────────────────────────────────────────────────────────────────
#   Load libraries
# ──────────────────────────────────────────────────────────────────────────────
library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# ──────────────────────────────────────────────────────────────────────────────
#   Functions
# ──────────────────────────────────────────────────────────────────────────────

# parseGENCODEgtf: Read a GENCODE GTF and extract attributes into columns
#
# Args:
#   gtf_path: path or URL to a GTF file (can be gzipped)
#   feature_types: optional character vector of feature_type values to keep (e.g., c("gene","transcript"))
#
# Returns:
#   A tibble with the standard GTF columns plus one column per attribute
parseGENCODEgtf <- function(gtf_path = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz", feature_types = c("gene","transcript")) {
  # 1. Read in the GTF
  message("1/6 ▶ Reading GTF...")
  raw_gtf <- read_tsv(
    gtf_path,
    comment = "##",
    col_names = c(
      "seqnames", "source", "feature_type",
      "start", "end", "score", "strand",
      "phase", "attributes"
    ),
    show_col_types = FALSE
  )

  # 2. Subset by feature_type if requested
  if (!is.null(feature_types)) {
    message("2/6 ▶ Subsetting by feature type...")
    raw_gtf <- raw_gtf %>%
      filter(feature_type %in% feature_types)
  }

  # 3. Define all possible attribute keys
  message("3/6 ▶ Identifying possible attribute keys...")
  # Hard coded
  # fields <- c(
  #   "gene_id", "gene_type", "gene_status", "gene_name",
  #   "transcript_id", "transcript_type", "transcript_status", "transcript_name",
  #   "exon_number", "exon_id", "level", "tag", "ccdsid",
  #   "havana_gene", "havana_transcript", "protein_id", "ont",
  #   "transcript_support_level", "remap_status", "remap_original_id",
  #   "remap_original_location", "remap_num_mappings", "remap_target_status",
  #   "remap_substituted_missing_target", "hgnc_id", "mgi_id"
  # )

  # Dynamic
  fields <- raw_gtf %>%
    separate_rows(attributes, sep = ";") %>%
    mutate(attributes = str_squish(attributes)) %>%
    filter(attributes != "") %>%
    mutate(key = word(attributes, 1)) %>%
    distinct(key) %>%
    pull(key)


  # 4. Compute max occurrence per attribute to identify single vs multi
  max_counts <- sapply(fields, function(fld) {
    counts_per_row <- str_count(
      raw_gtf$attributes,
      regex(paste0("\\b", fld, "\\b"))
    )
    max(counts_per_row, na.rm = TRUE)
  })

  # 5. Split into single- and multi- occurrence keys
  single_fields <- names(max_counts)[max_counts == 1]
  multi_fields  <- names(max_counts)[max_counts > 1]

  # 6. Extract single-occurrence attributes (NA if missing)
  message("4/6 ▶ Extracting single-occurrence keys...")
  for (f in single_fields) {
    pat <- paste0(f, ' "([^\"]+)"')
    raw_gtf[[f]] <- str_match(raw_gtf$attributes, pat)[,2]
  }

  # 7. Extract multi-occurrence attributes, collapse with '|'
  message("5/6 ▶ Extracting and collapsing multi-occurrence keys...")
  for (f in multi_fields) {
    pat <- paste0(f, ' "([^\"]+)"')
    raw_gtf[[f]] <- vapply(raw_gtf$attributes, function(attr) {
      vals <- str_match_all(attr, pat)[[1]][,2]
      if (length(vals) == 0) {
        NA_character_
      } else {
        paste(vals, collapse = "|")
      }
    }, FUN.VALUE = character(1))
  }

  # 8. Return the parsed tibble
  message("6/6 ▶ Returning formatted tibble...")
  raw_gtf
}


# annotateGTF: Add tag flags and external ID mappings to a parsed GTF tibble
#
# Args:
#   parsed_gtf: tibble returned by parseGENCODEgtf(), containing GTF columns and attribute fields
#   flags: character vector of tag names to flag in the “tag” column; defaults to
#          c("GENCODE_Primary", "Ensembl_canonical", "MANE_Select", "MANE_Plus_Clinical")
#   metadata_urls: named list of URLs for external metadata tables with elements
#                  - entrez: URL to Gencode EntrezGene metadata (tx → entrez_id)
#                  - hgnc:   URL to Gencode HGNC metadata (tx → HGNC_symbol, HGNC_id)
#                  - refseq: URL to Gencode RefSeq metadata (tx → RefSeq_TxID, RefSeq_PxID)
#                Defaults are set to Gencode release 48 locations.
#
# Returns:
#   A tibble that extends parsed_gtf with:
#     • one column per flag (0/1) for each tag in `flags`
#     • gene_id_base and transcript_id_base (version-stripped IDs)
#     • unique_id (row number)
#     • joined metadata columns: entrez_id, HGNC_symbol, HGNC_id, RefSeq_TxID, RefSeq_PxID
#     • TSS calculated from start/end and strand
annotateGTF <- function(
    parsed_gtf,
    flags = c(
      "GENCODE_Primary",
      "Ensembl_canonical",
      "MANE_Select",
      "MANE_Plus_Clinical"
    ),
    metadata_urls = list(
      entrez = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.metadata.EntrezGene.gz",
      hgnc   = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.metadata.HGNC.gz",
      refseq = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.metadata.RefSeq.gz"
    )
) {

  df <- parsed_gtf

  # 1) Add requested tag flags
  message("1/5 ▶ Adding tag flag fields...")
  if (length(flags) > 0) {
    for (f in flags) {
      df <- df %>%
        mutate(
          !!f := if_else(!is.na(tag) & str_detect(tag, fixed(f)), 1L, 0L)
        )
    }
  }

  # 2) Add base IDs and unique row identifier
  message("2/5 ▶ Adding base ensembl ids...")
  df <- df %>%
    mutate(
      gene_id_base       = sub("\\.\\d+$", "", gene_id),
      transcript_id_base = sub("\\.\\d+$", "", transcript_id),
      unique_id          = row_number()
    )

  # 3) Read & join metadata tables
  message("3/5 ▶ importing alternate id types from GENCODE metadata...")
  entrez <- read_tsv(metadata_urls$entrez,
                     col_names = c("transcript_id","entrez_id"),
                     show_col_types = FALSE)
  hgnc   <- read_tsv(metadata_urls$hgnc,
                     col_names = c("transcript_id","HGNC_symbol","HGNC_id"),
                     show_col_types = FALSE)
  refseq <- read_tsv(metadata_urls$refseq,
                     col_names = c("transcript_id","RefSeq_TxID","RefSeq_PxID"),
                     show_col_types = FALSE)

  message("4/5 ▶ Adding alternate id types...")
  df %>%
    left_join(entrez, by = "transcript_id", relationship = "many-to-many") %>%
    left_join(hgnc,   by = "transcript_id", relationship = "many-to-many") %>%
    left_join(refseq, by = "transcript_id", relationship = "many-to-many") %>%
    distinct()

  message("5/5 ▶ Calculating TSS...")
  df %>%
    mutate(
      TSS = if_else(strand == "+", start, end))

  }


# compute_match_stats: For each ID column, what fraction of `user_list` is found?
#
# Args:
#   user_list: character vector of IDs to map
#   map_df:    tibble containing ID‐columns to search
#   idCols:    character vector of column names to test (defaults to your standard set)
#
# Returns:
#   A data.frame with one row per id_type and columns:
#     • matched    — number of hits
#     • total      — length(user_list)
#     • pct_matched — matched/total * 100
#     • match_str   — “X of Y” label
compute_match_stats <- function(user_list,
                                map_df,
                                idCols = NULL) {
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


# map_and_report: Choose the best‐matching ID type for a user list and return minimal coordinates
#
# Args:
#   user_list:  character vector of gene or transcript IDs supplied by the user
#   map_df:     tibble produced by annotateGTF(), containing parsed GTF columns, flags, and external IDs
#   idCols:     optional character vector of column names in map_df to consider as possible ID matches;
#               if NULL, defaults to the standard set (gene_name, gene_id, transcript_id, etc.)
#   threshold:  numeric between 0 and 1 giving the minimum fraction of user_list that must match
#               (default 0.90) in the best‐hit column or the function returns NULL
#   gene_strategy: for gene‐level IDs, one of
#                 "gene", "flagged_transcript", "most_frequent_TSS",
#                 "unique_TSS", or "all_TSS"
#   flag_choice:   which flag to use when gene_strategy = "flagged_transcript"
#
# Returns:
#   A tibble with columns unique_id, seqnames, TSS, strand for the selected features,
#   or NULL if no ID type maps above the threshold.
map_and_report <- function(
    user_list,
    map_df,
    idCols        = NULL,
    threshold     = 0.90,
    gene_strategy = c(
      "gene",
      "flagged_transcript",
      "most_frequent_TSS",
      "unique_TSS",
      "all_TSS"
    ),
    flag_choice   = "GENCODE_Primary"
) {
  gene_strategy <- match.arg(gene_strategy)
  stats <- compute_match_stats(user_list, map_df, idCols)
  best  <- stats[which.max(stats$pct_matched), ]

  if (best$pct_matched / 100 < threshold) {
    message("Unable to map IDs ≥ ", threshold*100, "% in any field; please supply IDs in a supported format.")
    return(NULL)
  }
  message(sprintf(
    "Identified input ids as %s (%d/%d; %.1f%%).",
    best$id_type, best$matched, best$total, best$pct_matched
  ))

  # keep only rows matching user_list
  df <- map_df %>%
    filter(.data[[best$id_type]] %in% user_list)

  # classify level
  gene_keys       <- c("gene_id","gene_name","entrez_id","HGNC_symbol","HGNC_id")
  transcript_keys <- c("transcript_id","transcript_name","RefSeq_TxID")
  id_level <- if (best$id_type %in% gene_keys) "gene" else "transcript"

  if (id_level == "gene") {
    # apply gene_strategy...
    df_gene <- df %>% filter(feature_type %in% c("gene","transcript"))
    minimal <- switch(
      gene_strategy,
      gene = {
        df_gene %>% filter(feature_type == "gene")
      },
      flagged_transcript = {
        df_gene %>%
          filter(feature_type == "transcript", !!sym(flag_choice) == 1L)
      },
      most_frequent_TSS = {
        df_gene %>%
          filter(feature_type == "transcript") %>%
          mutate(TSS = if_else(strand == "+", start, end)) %>%
          count(gene_id, TSS, name = "n") %>%
          slice_max(n, with_ties = FALSE) %>%
          select(-n) %>%
          inner_join(df_gene %>% mutate(TSS = if_else(strand == "+", start, end)),
                     by = c("gene_id","TSS"))
      },
      unique_TSS = {
        df_gene %>%
          filter(feature_type == "transcript") %>%
          mutate(TSS = if_else(strand == "+", start, end)) %>%
          distinct(gene_id, TSS, .keep_all = TRUE)
      },
      all_TSS = {
        df_gene %>% filter(feature_type == "transcript")
      }
    )
  } else {
    message("Transcript-level IDs detected; returning all mapped transcripts (ignoring gene_strategy).")
    minimal <- df %>% filter(feature_type == "transcript")
  }

  # return only the promoter‐relevant columns, with unified TSS
  minimal %>%
    transmute(
      unique_id,
      seqnames,
      TSS    = if_else(strand == "+", start, end),
      strand
    )
}


# ──────────────────────────────────────────────────────────────────────────────
#   Code
# ──────────────────────────────────────────────────────────────────────────────
# Parse gencode gtf
gtf_path <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz"
gtf <- parseGENCODEgtf(feature_types = c("gene","transcript"))

# Add priority flags and alternate id types
gtf_map <- annotateGTF(gtf)


# Based on provided input list, map to an id type and return minimal table to use to identify promoter coordinates for matching genes/transcripts
user_list <- c(
  "BRCA1", "TP53", "ENST00000367770", "FOO", NA,
  "ENST00000832824.1", "ENST00000832825.1", "ENST00000832826.1",
  "ENST00000832827.1", "ENST00000832828.1", "ENST00000832829.1",
  "ENST00000832830.1", "ENST00000832837.1", "ENST00000832836.1",
  "ENST00000832832.1", "ENST00000832833.1", "ENST00000832831.1",
  "ENST00000832834.1", "ENST00000832835.1", "ENST00000832839.1",
  "ENST00000832841.1", "ENST00000832840.1", "ENST00000832838.1",
  "ENST00000832842.1", "ENST00000456328.3", "ENST00000456328.3",
  "ENST00000832845.1", "ENST00000832843.1", "ENST00000832844.1",
  "ENST00000832847.1", "ENST00000832846.1", "ENST00000832848.1",
  "ENST00000832849.1", "ENST00000832823.1", NA,
  "ENST00000450305.2", NA, "ENST00000831746.1", "ENST00000831201.1",
  "ENST00000831165.1", "ENST00000831154.1", "ENST00000831157.1",
  "ENST00000831145.1", "ENST00000831146.1", "ENST00000831177.1",
  "ENST00000831172.1", "ENST00000831158.1", "ENST00000831206.1",
  "ENST00000831205.1", "ENST00000831204.1", "ENST00000831185.1",
  "ENST00000831173.1", "ENST00000831170.1", "ENST00000831189.1"
)


user_list <- gtf_map$gene_id[1:1000]

Matched <- map_and_report(user_list,gtf_map,gene_strategy = "most_frequent_TSS")

