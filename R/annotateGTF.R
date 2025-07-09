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
      "GENCODE_Primary"
      ,"Ensembl_canonical"
      ,"MANE_Select"
      # ,"MANE_Plus_Clinical" # Probably not really a useful flag to include - it is very limited and doesn't include MANE_Select
    ),
    metadata_urls = list(
      entrez = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.metadata.EntrezGene.gz",
      hgnc   = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.metadata.HGNC.gz",
      refseq = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.metadata.RefSeq.gz"
    )
) {

  df <- parsed_gtf

  # # 1) Add requested tag flags
  # message("1/5 ▶ Adding tag flag fields...")
  # if (length(flags) > 0) {
  #   for (f in flags) {
  #     df <- df %>%
  #       mutate(
  #         !!f := if_else(!is.na(tag) & str_detect(tag, fixed(f)), 1L, 0L)
  #       )
  #   }
  # }

  # 1) Add requested tag flags (with “flag_” prefix)
  message("1/6 ▶ Adding tag flag fields...")
  if (length(flags) > 0) {
    for (f in flags) {
      flag_col <- paste0("flag_", f)
      df <- df %>%
        mutate(
          !!flag_col := as.integer(
            coalesce(
              str_detect(tag, fixed(f)),
              FALSE
            )
          )
        )
    }
  }

  # 2) Add base IDs
  message("2/6 ▶ Adding base ensembl ids...")
  df <- df %>%
    mutate(
      gene_id_base       = sub("\\.\\d+$", "", gene_id),
      transcript_id_base = sub("\\.\\d+$", "", transcript_id)
    )

  # 3) Read & join metadata tables
  message("3/6 ▶ importing alternate id types from GENCODE metadata...")
  entrez <- read_tsv(metadata_urls$entrez,
                     col_names = c("transcript_id","entrez_id"),
                     show_col_types = FALSE)
  hgnc   <- read_tsv(metadata_urls$hgnc,
                     col_names = c("transcript_id","HGNC_symbol","HGNC_id"),
                     show_col_types = FALSE)
  refseq <- read_tsv(metadata_urls$refseq,
                     col_names = c("transcript_id","RefSeq_TxID","RefSeq_PxID"),
                     show_col_types = FALSE)

  message("4/6 ▶ Adding alternate id types...")
  df %>%
    left_join(entrez, by = "transcript_id", relationship = "many-to-many") %>%
    left_join(hgnc,   by = "transcript_id", relationship = "many-to-many") %>%
    left_join(refseq, by = "transcript_id", relationship = "many-to-many") %>%
    distinct()

  message("5/6 ▶ Adding unique_id (from row name)...")
  df <- df %>%
    mutate(unique_id = row_number())

  message("6/6 ▶ Calculating TSS...")
  df %>%
    mutate(
      TSS = if_else(strand == "+", start, end))

}
