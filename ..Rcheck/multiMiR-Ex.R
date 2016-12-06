pkgname <- "multiMiR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('multiMiR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("get.multimir")
### * get.multimir

flush(stderr()); flush(stdout())

### Name: get.multimir
### Title: Get microRNA-target Interactions from the multiMiR Package
### Aliases: get.multimir
### Keywords: utilities database

### ** Examples

  ## search 'hsa-miR-18a-3p' in validated interactions in human
  example1 <- get.multimir(mirna='hsa-miR-18a-3p', summary=TRUE)
  names(example1)
  ## target genes that are validated by Luciferase assay
  example1$validated[grep("Luciferase", example1$validated[,"experiment"]),]
  example1$summary[example1$summary[,"target_symbol"] == "KRAS",]

  ## search 'cisplatin' in disease and drug tables in human
  example2 <- get.multimir(disease.drug='cisplatin', table='disease.drug')
  nrow(example2$disease.drug)
  head(example2$disease.drug)



cleanEx()
nameEx("list.multimir")
### * list.multimir

flush(stderr()); flush(stdout())

### Name: list.multimir
### Title: List microRNAs, Genes, Drugs Or Diseases in the multiMiR Package
### Aliases: list.multimir
### Keywords: utilities database

### ** Examples

  miRNAs <- list.multimir("mirna")
  genes <- list.multimir("gene")
  drugs <- list.multimir("drug")
  diseases <- list.multimir("disease")



cleanEx()
nameEx("multimir_dbInfo")
### * multimir_dbInfo

flush(stderr()); flush(stdout())

### Name: multimir_dbInfo
### Title: Collect Information About the Web Server And Database of the
###   multiMiR Package
### Aliases: getOption("multimir.url") multimir_dbCount multimir_dbInfo
###   multimir_dbSchema multimir_dbTables
### Keywords: utilities database

### ** Examples

  getOption("multimir.url")

  db.count <- multimir_dbCount()

  db.info <- multimir_dbInfo()

  multimir_dbSchema()

  db.tables <- multimir_dbTables()



cleanEx()
nameEx("search.multimir")
### * search.multimir

flush(stderr()); flush(stdout())

### Name: search.multimir
### Title: Search the multiMiR Database Given a MySQL Query
### Aliases: search.multimir
### Keywords: utilities database

### ** Examples

  ## show all tables in the multiMiR database
  tables <- search.multimir(query="show tables")

  ## show the structure of table diana_microt
  microt <- search.multimir(query="describe diana_microt")

  ## search for validated target genes of hsa-miR-18a-3p in miRecords
  result <- search.multimir(query="select m.mature_mirna_acc, m.mature_mirna_id, t.target_symbol, t.target_entrez, t.target_ensembl, i.experiment, i.support_type, i.pubmed_id from mirna AS m INNER JOIN mirecords AS i INNER JOIN target AS t ON (m.mature_mirna_uid=i.mature_mirna_uid and i.target_uid=t.target_uid) WHERE m.mature_mirna_id='hsa-miR-18a-3p'")



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
