################################################################################
# Script for testing during active development
################################################################################

# Load packages
library(devtools)
#library(roxygen2)
#library(Rd2roxygen)
library(dplyr)
library(tidyr)
library(purrr)


########################################
# Build/Check/Install/Load
########################################

setwd("~/Projects/KechrisLab/multiMiR")
devtools::document()
devtools::check()
devtools::install("../multiMiR", build_vignettes = TRUE)
library(multiMiR)

########################################
# Convert to roxygen2 documentation
########################################

# formatR::usage(Rd2roxygen)
# setwd("~/Projects/KechrisLab/")
# Rd2roxygen("multiMiR")

########################################
# Load old package for comparison/testing
# all function names use multimirold instead of multimir
########################################

install_github("mmulvahill/multiMiRold", ref = "refactor")
library(multiMiRold)

########################################
# Test/dev code
########################################

table <- "miranda"
mirna <- "mirna"
target <- "target"
org <- "rno"
mirna.table  <- "mirna.table"
target.table <- "target.table"
predicted.cutoff <- NULL
predicted.cutoff.type <- "n"
predicted.site <- "conserved"


# Check one single validated table
qry_tarbase <- get.multimir(table = "tarbase", mirna = "hsa-miR-199a-3p")
# Check all validated tables
qry_valid   <- get.multimir(table = "validated", mirna = "hsa-miR-199a-3p")

# Check one single disease.drug tables
qry_disease <- get.multimir(table = "mir2disease", mirna = "hsa-miR-199a-3p",
                            target = "TP53", org = "rno", 
                            disease.drug = "cisplatin")
# Check all disease.drug tables
qry_disease <- get.multimir(table = "disease.drug", mirna = "hsa-miR-199a-3p",
                            #target = "TP53",
                            org = "rno", disease.drug = "cisplatin")


##############################
# All possible options to get.multimir()
##############################

# TEST 1: Are organism names proper cleaned?
org <- c("hsa", "human", "homo sapiens", "mmu", "mouse", 
         "mus musculus", "rno", "rat", "rattus norvegicus")
# TEST 2: Is disease.drug case-insensitive? 
disease.drug = c(NULL, "bladder cancer", "cisplatin", "CiSplatin",
                 "blaDDER Cancer"),
# TEST 3: Do individual table queries combined equal combined categories?
table = c("validated", "predicted", "disease.drug", "all", "mirecords",
          "mirtarbase", "tarbase", "diana_microt", "elmmo", "microcosm",
          "miranda", "mirdb", "pictar", "pita", "targetscan", "mir2disease",
          "pharmaco_mir", "phenomir"),

get.multimir_opts <- 
    expand.grid(org = c("hsa", "mmu","rno"),
                # vector of mature miRNA accession #s or mature miRNA IDs, or
                # combo just an example
                mirna = c(NULL, "MIMAT0000072", "hsa-miR-199a-3p",
                          "MIMAT0000065", "hsa-miR-30a-5p"),
                # target -- vector of gene symbols, Entrez gene IDs, Ensembl
                # gene IDs, or combo -- just an example
                target = c(NULL, "TP53", "KRAS", 578, 3845, "ENSG00000171791",
                           "TP53", 3845, "ENSG00000171791"),
                # disease or drug -- vector of diseases and/or drugs
                # just an example, should be case-insensitive
                disease.drug = c(NULL, "bladder cancer", "cisplatin"),
                # Table type or table name -- all possible options, if table
                # name only 1 allowed
                table = c("validated", "predicted", "disease.drug"),
                # Prediction score cutoffs -- 0:100 if p=percent, or > 10000 if
                # n=count
                predicted.cutoff = c(NULL, seq(0, 100, by = 10), 
                                     seq(10000, 500000, by = 10000)),
                predicted.cutoff.type = c(NULL, "p", "n"), # --all possible options
                # type of predicted target site to search -- only works for
                # miranda, pita, and targetscan
                predicted.site = c("conserved", "nonconserved", "all")
            )

# non-query options
summary <- c(TRUE, FALSE)
add.link <- c(TRUE, FALSE)
url <- NULL





