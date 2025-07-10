library(readr)    # for read_tsv()
library(dplyr)    # for filter(), mutate(), distinct(), pull()
library(tidyr)    # for separate_rows()
library(stringr)  # for str_squish(), str_count(), str_match(), str_match_all(), word()


# parseGENCODEgtf: Read a GENCODE GTF and extract attributes into columns
#
# Args:
#   gtf_path: path or URL to a GTF file (can be gzipped)
#   feature_types: optional character vector of feature_type values to keep (e.g., c("gene","transcript"))
#
# Returns:
#   A tibble with the standard GTF columns plus one column per attribute
parseGENCODEgtf <- function(gtf_path = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz",
                            feature_types = c("gene","transcript")) {
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
