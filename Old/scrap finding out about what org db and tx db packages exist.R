
# get the Bioconductor repo URLs
repos <- BiocManager::repositories()

# pull the available package list from CRAN+Bioconductor
pkgs  <- available.packages(repos = repos)

# filter to all org.*.eg.db packages
orgdbs <- grep("^org\\.[A-Za-z]+\\.eg\\.db$", rownames(pkgs), value = TRUE)


orgdbASCII <- grep("^org\\.[A-Za-z]+\\.[A-Za-z]+\\.db$", rownames(pkgs), value = TRUE)
orgdbAll <- grep("^org\\.[^.]+\\.[^.]+\\.db$", rownames(pkgs), value=TRUE)

orgdbASCII
orgdbAll

setdiff(orgdbAll,orgdbASCII)

> orgdbAll
# [1] "org.Ag.eg.db"      "org.At.tair.db"    "org.Bt.eg.db"      "org.Ce.eg.db"      "org.Cf.eg.db"      "org.Dm.eg.db"      "org.Dr.eg.db"
# [8] "org.EcK12.eg.db"   "org.EcSakai.eg.db" "org.Gg.eg.db"      "org.Hs.eg.db"      "org.Mm.eg.db"      "org.Mmu.eg.db"     "org.Pf.plasmo.db"
# [15] "org.Pt.eg.db"      "org.Rn.eg.db"      "org.Sc.sgd.db"     "org.Ss.eg.db"      "org.Xl.eg.db"

# Only ones not using entrez geneid are Arabidopsis thaliana, Plasmodium falciparum, Saccharomyces cerevisiae (yeast)


BiocManager::install("org.At.tair.db")
library(org.At.tair.db)



# see all keytypes available
keytypes(org.At.tair.db)

tair_ids <- keys(org.At.tair.db, keytype="TAIR")
ENTREZ_ids <- keys(org.At.tair.db, keytype = "ENTREZID")






# 3. Filter to any package whose name starts with “TxDb.”
txdb_pkgs <- grep("^TxDb\\.", rownames(pkgs), value = TRUE)

# 4. Inspect
length(txdb_pkgs)
# e.g. [1] 14
txdb_pkgs
# [1] "TxDb.Athaliana.BioMart.plantsmart25"  "TxDb.Celegans.UCSC.ce11.refGene"
#      "TxDb.Dmelanogaster.UCSC.dm6.ensGene" ... etc.





library(TxDb.Hsapiens.UCSC.hg38.knownGene)



BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

# 1a) List all the “keytypes” you can query by
kt <- keytypes(org.Hs.eg.db)
print(kt)
# e.g. "ENTREZID", "SYMBOL", "REFSEQ", "ENSEMBL", "UNIPROT", "GO", …

# 1b) List all the “columns” you can pull out
cols <- AnnotationDbi::columns(org.Hs.eg.db)
print(cols)


# e.g. matches keytypes, plus some extras like "GENENAME", "PATH", "PMID", …

# 1c) Show how many entries map in each
sapply(kt, function(k) length(keys(org.Hs.eg.db, keytype = k)))


# 1) Grab the first 20 Entrez Gene IDs
entrez_keys <- head(
  keys(org.Hs.eg.db, keytype="ENTREZID"),
  20
)


df <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = entrez_keys,
  keytype = "ENTREZID",
  columns = cols
)

# 3) Inspect
head(df)

# 1. Get the combined CRAN+Bioconductor repo URLs
repos <- BiocManager::repositories()

# 2. Pull the current package list
pkgs <- available.packages(repos = repos)

# 3. Filter to any package whose name starts with “TxDb.”
txdb_pkgs <- grep("^TxDb\\.", rownames(pkgs), value = TRUE)

# 4. Inspect
length(txdb_pkgs)
# e.g. [1] 14
txdb_pkgs
# [1] "TxDb.Athaliana.BioMart.plantsmart25"  "TxDb.Celegans.UCSC.ce11.refGene"
#      "TxDb.Dmelanogaster.UCSC.dm6.ensGene" ... etc.


all_bioc <- BiocManager::available()
txdb_pkgs <- grep("^TxDb\\.", all_bioc, value = TRUE)



BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
