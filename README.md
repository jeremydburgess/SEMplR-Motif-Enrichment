# GENCODE Promoter-Extraction Workflow Primer

This primer describes the end-to-end pipeline for extracting promoter regions from GENCODE annotations, mapping user-supplied gene/transcript IDs, and preparing foreground TSS windows for motif enrichment.

---

## 1. parseGENCODEgtf.R → `parseGENCODEgtf()`

**Purpose:**

- Read a (gzipped) GENCODE GTF, subset to specified feature types (`gene`, `transcript`), and split the final `attributes` column into one tidy column per key.
- Collapse multi-valued keys (e.g. `tag`, `ont`) with `|`.

**Inputs:**

- `gtf_path`: file path or URL to `.gtf(.gz)`
- `feature_types`: character vector to keep (default `c("gene","transcript")`)

**Output:**

- Tibble with the 8 standard GTF columns plus one column per attribute key.

---

## 2. annotateGTF.R → `annotateGTF()`

**Purpose:**

- Extend the parsed GTF tibble with:
  - Boolean `flag_<TAG>` columns for selected tags (default: `GENCODE_Primary`, `Ensembl_canonical`, `MANE_Select`, `MANE_Plus_Clinical`).
  - `gene_id_base` and `transcript_id_base` (version-stripped IDs).
  - `unique_id` (row number).
  - External IDs via left joins to GENCODE Entrez, HGNC, and RefSeq metadata.
  - A unified `TSS` column computed from `start`/`end` and `strand`.

**Inputs:**

- `parsed_gtf`: output of `parseGENCODEgtf()`
- `flags`: vector of tags to flag (defaults above)

**Output:**

- Tibble with all parsed columns plus the new annotations and `TSS`.

---

## 3. computeMatchStats.R → `compute_match_stats()`

**Purpose:**

- For each candidate ID column, calculate how many of the user’s IDs match and the percentage mapped.

**Inputs:**

- `user_list`: character vector of IDs to map
- `map_df`: annotated GTF tibble
- `idCols`: optional override of columns to test

**Output:**

- Data frame with rows per `id_type`, and columns `matched`, `total`, `pct_matched`, `match_str`.

---

## 4. mapUserList.R → `mapUserList()`

**Purpose:**

- Run `compute_match_stats()`, pick the best-matching column, determine `id_level` (`gene` vs `transcript`), and subset the annotated GTF to only those matching rows.
- Add `mappedInput_id` (the original matched ID) and ensure `TSS` exists.

**Inputs:**

- `user_list`, `map_df`, optional `idCols`, `threshold` (default 0.9)

**Output:**

- List with elements:
  - `match_stats`, `best_id_type`, `id_level`, `mapped_df` (the subset tibble).

---

## 5. assembleTSSforeground.R → `assembleTSSforeground()`

**Purpose:**

- From the `mapUserList()` output, build a foreground TSS table according to a `gene_filter` strategy:
  - `gene`: keep gene entries.
  - `flag:<TAG>`: keep transcripts with `flag_<TAG> == 1`.
  - `modalTSS`: one TSS per gene at the most-frequent coordinate.
  - `uniqueTSS`: one row per distinct TSS per gene.
  - `allTSS`: all transcript TSSs.

**Inputs:**

- `mapped_data` (list from `mapUserList()`)
- `id_level = c("auto","gene","transcript")`
- `gene_filter` (see above)

**Output:**

- Tibble of `unique_id`, `seqnames`, `TSS`, `strand` (and any flag columns) ready for promoter region calculation.

---

## 6. definePromoterRegions.R → `definePromoterRegions()`

**Purpose:**

- Convert a TSS tibble into a **GRanges** of promoter windows (e.g. -1 kb/+100 bp), preserving metadata (e.g. `gene_id`, `unique_id`).

**Inputs:**

- Tibble with `seqnames`, `strand`, `TSS`, plus extra columns.
- `region = c(upstream, downstream)` offsets in bp.

**Output:**

- `GRanges` object with ranges `[TSS–up, TSS+down]` (strand-aware) and metadata columns.

---

## 7. RunPipeline.R (Example Script)

Demonstrates end-to-end usage:

```r
source("parseGENCODEgtf.R")
source("annotateGTF.R")
source("computeMatchStats.R")
source("mapUserList.R")
source("assembleTSSforeground.R")
source("definePromoterRegions.R")

gtf_raw  <- parseGENCODEgtf(gtf_path)
gtf_map  <- annotateGTF(gtf_raw)
# Randomly select 500 gene_ids as the user-supplied list
set.seed(123)  # for reproducibility
user_list <- gtf_map$gene_id[sample.int(nrow(gtf_map), 500, replace = FALSE)]

#  Run mapping and assemble promoter regions
mapped   <- mapUserList(user_list, gtf_map)
fg_tss   <- assembleTSSforeground(mapped)
fg_prom  <- definePromoterRegions(fg_tss, region = c(300,50))

```
