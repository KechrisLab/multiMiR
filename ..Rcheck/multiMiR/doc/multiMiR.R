### R code from vignette source 'multiMiR.Rnw'

###################################################
### code chunk number 1: annotate
###################################################
options(width=100)


###################################################
### code chunk number 2: multimir_dbTables
###################################################
library(multiMiR)
db.tables = multimir_dbTables()
db.tables


###################################################
### code chunk number 3: multiMiR.Rnw:50-51 (eval = FALSE)
###################################################
## multimir_dbSchema()


###################################################
### code chunk number 4: multimir_dbInfo
###################################################
db.info = multimir_dbInfo()
db.info


###################################################
### code chunk number 5: multimir_dbCount
###################################################
db.count = multimir_dbCount()
db.count
apply(db.count[,-1], 2, sum)


###################################################
### code chunk number 6: list.multimir
###################################################
miRNAs = list.multimir("mirna")
dim(miRNAs)
head(miRNAs)
genes = list.multimir("gene")
dim(genes)
head(genes)
drugs = list.multimir("drug")
dim(drugs)
head(drugs)
diseases = list.multimir("disease")
dim(diseases)
head(diseases)


###################################################
### code chunk number 7: Example1
###################################################
# The default is to search validated interactions in human
example1 = get.multimir(mirna='hsa-miR-18a-3p', summary=TRUE)
names(example1)
# Detailed information of the validated miRNA-target interaction
head(example1$validated)
# Which interactions are supported by Luciferase assay?
example1$validated[grep("Luciferase", example1$validated[,"experiment"]),]
example1$summary[example1$summary[,"target_symbol"] == "KRAS",]


###################################################
### code chunk number 8: Example2
###################################################
example2 = get.multimir(disease.drug='cisplatin', table='disease.drug')
names(example2)
nrow(example2$disease.drug)
head(example2$disease.drug)


###################################################
### code chunk number 9: Example3_part1
###################################################
example3 = get.multimir(org="mmu", target="Gnb1", table="predicted", summary=TRUE,
predicted.cutoff=35, predicted.cutoff.type="p", predicted.site="all")
names(example3)
head(example3$predicted)
head(example3$summary)


###################################################
### code chunk number 10: Example3_part2
###################################################
apply(example3$summary[,6:13], 2, function(x) sum(x>0))


###################################################
### code chunk number 11: Example4_part1
###################################################
example4 = get.multimir(org='hsa', target=c('AKT2','CERS6','S1PR3','SULF2'), table='predicted',
summary=TRUE, predicted.cutoff.type='n', predicted.cutoff=500000)


###################################################
### code chunk number 12: Example4_part2
###################################################
example4.counts = addmargins(table(example4$summary[,2:3]))
example4.counts = example4.counts[-nrow(example4.counts),]
example4.counts = example4.counts[order(example4.counts[,5], decreasing=TRUE),]
head(example4.counts)


###################################################
### code chunk number 13: multiMiR.Rnw:259-260 (eval = FALSE)
###################################################
## load(url("http://multimir.ucdenver.edu/bladder.rda"))


###################################################
### code chunk number 14: multiMiR.Rnw:267-270 (eval = FALSE)
###################################################
## # search all tables & top 10% predictions
## example5 = get.multimir(org="hsa", mirna=DE.miRNA.up, target=DE.entrez.dn, table="all",
## summary=TRUE, predicted.cutoff.type="p", predicted.cutoff=10)


###################################################
### code chunk number 15: multiMiR.Rnw:293-294 (eval = FALSE)
###################################################
## example5$validated


###################################################
### code chunk number 16: multiMiR.Rnw:315-316 (eval = FALSE)
###################################################
## example5$disease.drug[grep("bladder", example5$disease.drug$disease_drug, ignore.case=TRUE),]


###################################################
### code chunk number 17: multiMiR.Rnw:337-338 (eval = FALSE)
###################################################
## length(unique(example5$predicted$mature_mirna_id))


###################################################
### code chunk number 18: multiMiR.Rnw:343-344 (eval = FALSE)
###################################################
## length(unique(example5$predicted$target_entrez))


###################################################
### code chunk number 19: multiMiR.Rnw:349-352 (eval = FALSE)
###################################################
## unique.pairs = unique(data.frame(miRNA.ID=as.character(example5$predicted$mature_mirna_id),
## target.Entrez=as.character(example5$predicted$target_entrez)))
## nrow(unique.pairs)


###################################################
### code chunk number 20: multiMiR.Rnw:357-358 (eval = FALSE)
###################################################
## head(unique.pairs)


###################################################
### code chunk number 21: multiMiR.Rnw:374-375 (eval = FALSE)
###################################################
## example5.split = split(example5$predicted, example5$predicted$database)


###################################################
### code chunk number 22: Direct_query1
###################################################
direct1 = search.multimir(query="show tables")


###################################################
### code chunk number 23: Direct_query2
###################################################
direct2 = search.multimir(query="describe diana_microt")
direct2


###################################################
### code chunk number 24: Direct_query3
###################################################
direct3 = search.multimir(query="select m.mature_mirna_acc, m.mature_mirna_id, t.target_symbol, t.target_entrez, t.target_ensembl, i.experiment, i.support_type, i.pubmed_id from mirna AS m INNER JOIN mirecords AS i INNER JOIN target AS t ON (m.mature_mirna_uid=i.mature_mirna_uid and i.target_uid=t.target_uid) WHERE m.mature_mirna_id='hsa-miR-18a-3p'")
direct3


###################################################
### code chunk number 25: sessionInfo
###################################################
sessionInfo()
warnings()


