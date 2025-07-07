# # Use Biocmanager packages to generate a flat file to convert between ID types
#
# # Type	Example	Use Case
# # HGNC Symbol	TP53	Most human-readable; standard in papers
# # Ensembl Gene ID	ENSG00000141510	Consistent, transcript-agnostic
# # Ensembl Transcript ID	ENST00000269305	Transcript-specific (e.g. isoforms)
# # Entrez Gene ID	7157	Often used in legacy datasets/tools
# # RefSeq Gene ID	NM_000546 or NR_123456	Transcript-level from NCBI
# # UniProt ID	P04637 (for proteins)	Often in proteomics or TF databases
#
#
# # ID Mapping Script: Build a reusable cross-reference table for gene/transcript identifiers
#
# # BiocManager::install("clusterProfiler")
# # BiocManager::install("org.Hs.eg.db")
#
#
# # Load necessary libraries
# library(rtracklayer)       # For GTF import
# library(data.table)        # For fast reading of large tables
# library(dplyr)             # Data manipulation
# library(tidyr)             # Data cleaning
# library(clusterProfiler)   # For bitr conversion
# library(org.Hs.eg.db)      # OrgDb for human gene annotation
#
# # ---- Step 1: Import GENCODE GTF ----
gtf_path <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz"  # Change to your GTF path
gtf <- import(gtf_path)

# Extract gene-level info
genes <- gtf[gtf$type == "gene"]
gencode_genes <- as.data.frame(mcols(genes)) %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  distinct()

# # ---- Step 2: Download and parse NCBI gene info ----
# ncbi_url <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
# ncbi <- fread(ncbi_url)
#
# ncbi_subset <- ncbi[, .(GeneID, Symbol, Synonyms, dbXrefs)] %>%
#   rename(
#     Entrez = GeneID,
#     NCBI_Symbol = Symbol,
#     NCBI_Aliases = Synonyms,
#     NCBI_DBXrefs = dbXrefs
#   )
#
# # ---- Step 3: Use org.Hs.eg.db for SYMBOL ↔ ENSEMBL/ENTREZ/ALIAS mapping ----
# orgdb_map <- bitr(unique(gencode_genes$gene_name),
#                   fromType = "SYMBOL",
#                   toType = c("ENSEMBL", "ENTREZID", "ALIAS"),
#                   OrgDb = org.Hs.eg.db) %>%
#   rename(
#     HGNC_Symbol = SYMBOL,
#     EnsemblGene = ENSEMBL,
#     Entrez = ENTREZID,
#     Alias = ALIAS
#   )
#
# # ---- Step 4: UCSC kgAlias.txt.gz ----
# ucsc_alias_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/kgAlias.txt.gz"
# ucsc_alias <- fread(ucsc_alias_url, header = FALSE)
# colnames(ucsc_alias) <- c("UCSC_Transcript", "UCSC_Alias")
#
# # ---- Optional: Transcript ↔ Gene map from GENCODE ----
tx <- gtf[gtf$type == "transcript"]
tx_map <- as.data.frame(mcols(tx)) %>%
  dplyr::select(transcript_id, gene_id, gene_name, transcript_type) %>%
  distinct()
#
# # ---- Combine all ----
# # Start with GENCODE genes as base
# id_map <- gencode_genes %>%
#   left_join(orgdb_map, by = c("gene_name" = "HGNC_Symbol")) %>%
#   left_join(ncbi_subset, by = c("gene_name" = "NCBI_Symbol"))
#
# # Save full mapping
# write_tsv(id_map, "ID_mapping_master.tsv")
#
# # Optional: Save transcript-level map too
# write_tsv(tx_map, "GENCODE_transcript_map.tsv")




# Install if not already installed
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}

# Load the package
library(biomaRt)


# Use Ensembl's current human dataset (GRCh38 / hg38)
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")


# Define the attributes you want to retrieve
attributes <- c(
  "ensembl_transcript_id",
  "ensembl_gene_id",
  "external_gene_name",
  "transcript_biotype",
  "refseq_mrna",
  "chromosome_name",
  "transcript_start",
  "transcript_end",
  "strand"
)

# Get the data
transcripts <- getBM(attributes = attributes, mart = ensembl)




# # ---- Step 1: Import GENCODE GTF ----
gtf_path <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz"  # Change to your GTF path
gtf <- rtracklayer::import(gtf_path)
gtf_frame <- as.data.frame(gtf)

colnames(gtf_frame)

# # Extract gene-level info
# genes <- gtf[gtf$type == "gene"]
# gencode_genes <- as.data.frame(mcols(genes)) %>%
#   dplyr::select(gene_id, gene_name, gene_type, hgnc_id, havana_gene) %>%
#   distinct()


# # ---- Optional: Transcript ↔ Gene map from GENCODE ----
tx <- gtf[gtf$type == "transcript"]
tx_map <- as.data.frame(mcols(tx)) %>%
  dplyr::select(transcript_id, transcript_type, gene_id, gene_name, gene_type, havana_gene)


entrez_meta_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.metadata.EntrezGene.gz"
HGNC_meta_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.metadata.HGNC.gz"
RefSeq_meta_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.metadata.RefSeq.gz"

entrez_meta <- read_tsv(entrez_meta_url, col_names = c("transcript_id", "entrez_id"), show_col_types = FALSE)
HGNC_meta <- read_tsv(HGNC_meta_url, col_names = c("transcript_id", "HGNC_symbol", "HGNC_id"), show_col_types = FALSE)
RefSeq_meta <- read_tsv(RefSeq_meta_url, col_names = c("transcript_id", "RefSeq_TxId", "RefSeq_PxId"), show_col_types = FALSE)

library(dplyr)

tx_annotated <- tx_map %>%
  left_join(entrez_meta, by = "transcript_id") %>%
  left_join(HGNC_meta, by = "transcript_id") %>%
  left_join(RefSeq_meta, by = "transcript_id") %>%
  distinct()

gene_map <- tx_annotated %>%
  dplyr::select(gene_id, gene_name, gene_type, havana_gene, entrez_id, HGNC_symbol, HGNC_id) %>%
  distinct()
