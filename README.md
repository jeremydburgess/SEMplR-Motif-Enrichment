## SEMplR-Motif-Enrichment

**Purpose:**  
Provide a streamlined pipeline for taking an arbitrary list of differentially expressed IDs,
mapping them to a universal ID space, extracting promoter regions (with overlap‐handling,
background selection, etc.), and handing off to SEMplR’s core motif‐enrichment routines.

---

### Pipeline Overview

1. **buildMappingObject()**  
   Instantiate an `src_organism` object by loading the appropriate `OrgDb` and `TxDb` for the user-specified organism.

2. **mapForegroundIds()**  
   - Detect which user‐supplied ID type (e.g. SYMBOL, ENSEMBL, etc.) matches best.  
   - Map to Entrez IDs.  
   - Construct initial “foreground” and “background universe” vectors of Entrez IDs.

3. **poolFilter()**  
   Subset the universe (foreground + background) to a desired gene type (e.g. protein-coding) or transcript feature, based on user settings.

4. **getCoordinates()**  
   - Pull genomic ranges from the TxDb.  
   - Define promoters (e.g. ±2 kb from TSS, or user-defined).  
   - Handle overlaps (merge or trim).  
   - Reduce multiple transcripts per gene to a single range (if requested).  
   - Sample a background set of promoters matched by length/GC content.  
   - Return final `GRanges` objects:  
     - `fgRanges` (foreground promoters)  
     - `bgRanges` (background promoters)

5. **enrichmentSets()**  
   A wrapper that ties all of the above together and returns a list of `GRanges` ready for SEMplR’s enrichment engine.

---

### Exported Functions

Below is a quick list of the key user-facing functions in **R/**:

| Function               | Description                                                                                             |
|------------------------|---------------------------------------------------------------------------------------------------------|
| `buildMappingObject()` | Given `organism` (e.g. `"Homo sapiens"`) and optional database paths, loads OrgDb & TxDb into a single `src_organism`. |
| `mapForegroundIds()`   | Maps your DE IDs to EntrezIDs, picks best ID type, and builds initial foreground & background ID vectors. |
| `poolFilter()`         | Filters your universe to a specific gene/transcript category (e.g. `"protein_coding"`).                |
| `getCoordinates()`     | Extracts promoter `GRanges` for each gene, resolves overlaps, selects background matches.               |
| `enrichmentSets()`     | High-level wrapper that runs the full pipeline and returns a named list `list(fg = fgRanges, bg = bgRanges)`. |

---

### Usage Example

```r
library(SEMplR.Motif.Enrichment)

# 1) Create the mapping object for human
so <- buildMappingObject(
  organism = "Homo sapiens",
  txdb     = "TxDb.Hsapiens.UCSC.hg38.knownGene",
  orgdb    = "org.Hs.eg.db"
)

# 2) Map your DE results (e.g. gene symbols) to Entrez + background
mapped <- mapForegroundIds(
  so,
  foreground_ids = c("BRCA1", "TP53", "MYC"),
  id_type        = NULL    # auto‐detect
)

# 3) Filter to protein-coding genes only
filtered <- poolFilter(
  mapped,
  gene_type = "protein_coding"
)

# 4) Get promoter GRanges (±2 kb around TSS)
coords <- getCoordinates(
  so,
  filtered,
  promoter_upstream   = 2000,
  promoter_downstream = 2000,
  overlap_action      = "reduce",
  collapse_method     = "onePerGene",
  n_background        = 1000
)

# coords$fg  # foreground GRanges
# coords$bg  # background GRanges

# 5) Or just do it all in one go:
sets <- enrichmentSets(
  organism           = "Homo sapiens",
  foreground_ids     = my_de_gene_list,
  gene_type          = "protein_coding",
  promoter_upstream   = 2000,
  promoter_downstream = 2000,
  n_background        = 1000
)
```

---

### Function Details

#### `buildMappingObject(organism, txdb, orgdb, ...)`  
- **organism**: full species name, e.g. `"Mus musculus"`  
- **txdb** / **orgdb**: package names (or file paths) for TxDb and OrgDb  
- Loads both and returns an S4 object containing both DB connections.

#### `mapForegroundIds(so, foreground_ids, id_type = NULL, ...)`  
- **so**: the object from `buildMappingObject()`  
- **foreground_ids**: character vector of IDs from your DE analysis  
- **id_type**: if `NULL`, auto‐detect among keys in the OrgDb  
- Returns a list with at least:  
  - `fg_entrez` (mapped foreground Entrez IDs)  
  - `bg_entrez` (all Entrez in the universe)

#### `poolFilter(mapped, gene_type = NULL, transcript_type = NULL)`  
- **mapped**: output of `mapForegroundIds()`  
- **gene_type** / **transcript_type**: e.g. `"protein_coding"`, vector of types  
- Returns a filtered version of the mapped list, with `fg_entrez` and `bg_entrez` restricted.

#### `getCoordinates(so, mapped, promoter_upstream, promoter_downstream, overlap_action = c("reduce","trim"), collapse_method = c("onePerGene","keepAll"), n_background = 1000, ...)`  
- **promoter_upstream/downstream**: numeric in bp  
- **overlap_action**: how to handle overlapping ranges  
- **collapse_method**: how to reduce multiple transcripts per gene  
- **n_background**: number of background ranges to sample, matched on length/GC  
- Returns a list:  
  - `fgRanges`: `GRanges` for foreground promoters  
  - `bgRanges`: `GRanges` for matched background promoters

#### `enrichmentSets(...)`  
- A convenience wrapper. All of the above arguments can be passed directly to this single function to get your final `fg` and `bg` GRanges.

---

---

*(This module is under active development. If you bump into any edge cases or want to tweak the filtering or promoter definitions, please open an issue or PR in this repo.)*
