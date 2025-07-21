helper










library(AnnotationHub)
ah <- AnnotationHub()
ensdb_records <- query(ah, c("EnsDb","homo sapiens"))
mcols(ensdb_records)$genome





library(GenomeInfoDb)
genomeStyles()
listOrganisms()
genomeBuilds("Homo sapiens", style = c("UCSC", "ensemblId"))

listOrganisms()

genomeBuilds("mouse")
genomeBuilds(c("Mouse", "dog", "human"), style=c("Ensembl","UCSC"))

mapGenomeBuilds(c("canFam3", "GRCm38", "mm9"))
mapGenomeBuilds(c("canFam3", "GRCm38", "mm9"), style="Ensembl")


human <- genomeBuilds("Homo sapiens")
mapGenomeBuilds(human$ucscID,style = "Ensembl")[,c("ensemblID","ucscID")]


UCSC <- genomeBuilds("Homo sapiens", style = "UCSC")
Ensembl <- genomeBuilds("Homo sapiens", style = "Ensembl")


# find the full path to the CSV
csv_file <- system.file(
  "extdata", "dataFiles", "genomeMappingTbl.csv",
  package = "GenomeInfoDb"
)
# read it in
genomeMap <- read.csv(csv_file, stringsAsFactors = FALSE)



# look at the first few titles
head(ensdb_records$title)

humanEnsDb <- data.frame(ensdb_records$ah_id,ensdb_records$title)

edb <- ensdb_records[["AH119325"]]  # replace AH#### with the record ID


library(ensembldb)

# 1) What tables are available?
listTables(edb)
# e.g. "transcripts", "genes", "exons", "cds", "promoters", etc.

# 2) What columns does a particular table have?
columns(edb, table="transcripts")
# or for genes:
listColumns(edb, table="gene")
listColumns(edb, table="tx")


genes(edb)
transcripts(edb)
# 3) High‐level metadata
metadata(edb)


