################################################################################
# Script for testing during active development
################################################################################

# Load packages
library(devtools)
library(testthat)
#library(roxygen2)
#library(Rd2roxygen)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)


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

# install_github("mmulvahill/multiMiRold", ref = "refactor")
library(multiMiRold)

########################################
# Dev code
########################################
# table <- "miranda"
# mirna <- "mirna"
# target <- "target"
# org <- "rno"
# mirna.table  <- "mirna.table"
# target.table <- "target.table"
# predicted.cutoff <- NULL
# predicted.cutoff.type <- "n"
# predicted.site <- "conserved"

########################################
# Test code
########################################
# Check one single validated table
# qry_tarbase_old <- get.multimirold(table = "tarbase", mirna = "hsa-miR-199a-3p")



#qry_tarbase     <- get.multimir(table = "tarbase", mirna = "hsa-miR-199a-3p")



# Check all validated tables
# qry_valid_old   <- get.multimirold(table = "validated", mirna = "hsa-miR-199a-3p")
# qry_valid       <- get.multimir(table = "validated", mirna = "hsa-miR-199a-3p")
# 
# # Check one single disease.drug tables
# qry_disease_old <- get.multimirold(table = "mir2disease", 
#                                    mirna = "hsa-miR-199a-3p",
#                                    target = "TP53", org = "rno", 
#                                    disease.drug = "cisplatin")
# # Check one single disease.drug tables
# qry_disease     <- get.multimir(table = "mir2disease", mirna = "hsa-miR-199a-3p",
#                                 target = "TP53", org = "rno", 
#                                 disease.drug = "cisplatin")
# # Check all disease.drug tables
# qry_disease     <- get.multimir(table = "disease.drug", mirna = "hsa-miR-199a-3p",
#                                 #target = "TP53",
#                                 org = "rno", disease.drug = "cisplatin")


##############################
# All possible options to get.multimir()
##############################

# # TEST 2: Is disease.drug case-insensitive? 
# disease.drug <- c(NULL, "bladder cancer", "cisplatin", "CiSplatin",
#                   "blaDDER Cancer"),
# # TEST 3: Do individual table queries combined equal combined categories?
# table <- c("validated", "predicted", "disease.drug", "all", "mirecords",
#            "mirtarbase", "tarbase", "diana_microt", "elmmo", "microcosm",
#            "miranda", "mirdb", "pictar", "pita", "targetscan", "mir2disease",
#            "pharmaco_mir", "phenomir"),
# 
# args_get.multimir <- 
#     expand.grid(org = c("hsa", "mmu","rno"),
#                 # vector of mature miRNA accession #s or mature miRNA IDs, or
#                 # combo -- just an example
#                 mirna = c(NULL, "MIMAT0000072", "hsa-miR-199a-3p"),
#                 # target -- vector of gene symbols, Entrez gene IDs, Ensembl
#                 # gene IDs, or combo -- just an example
#                 target = c(NULL, "TP53", 578, "ENSG00000171791"),
#                 # disease or drug -- vector of diseases and/or drugs
#                 # just an example, should be case-insensitive
#                 disease.drug = c(NULL, "bladder cancer", "cisplatin"),
#                 # Table type or table name -- all possible options, if table
#                 # name only 1 allowed
#                 table = c("validated", "predicted", "disease.drug")
#             ) %>% tbl_df
# 
# # Prediction score cutoffs -- 0:100 if p=percent, or > 10000 if
# # n=count
# 
# pred <- 
#     data.frame(predicted.cutoff      = c(0, 50, 100, 101, 9999, 10000, 500000),
#                predicted.cutoff.type = c(rep("p", 4), rep("n", 3)), 
#                # type of predicted target site to search -- only works for
#                # miranda, pita, and targetscan
#                predicted.site = rep(c("conserved", "nonconserved", "all"), 
#                                     each = 7)) %>% tbl_df
# 
# # non-query options
# summary <- c(TRUE, FALSE)
# add.link <- c(TRUE, FALSE)
# url <- NULL
# 

########################################
# Test functions for new vs old pkg
########################################
options(stringsAsFactors = FALSE)

test_arg_df <- function(arg_df) {
    arg_df %>%
        by_row(to_arglist, .to = 'arglist') %>% 
        mutate(queries   = map(arglist, qryget)) %>%
        mutate(success   = map(queries, qryequal)) %>% 
        mutate(allsucess = map(success, ~ all(as.logical(.x))) %>% unlist)
}
to_arglist <- function(dfrow) {
    dfrow %>% flatten
}
qryget <- function(getargs) {
    newfn <- do.call(get.multimir, getargs)
    oldfn <- do.call(get.multimirold, getargs)[[1]]
    list(newfn, oldfn)
}
qryequal <- function(newold_list) {
    new <- newold_list[[1]]
    old <- newold_list[[2]]
    if (length(old) == length(new)) {
        rtn <- lapply(1:length(old), compare_qrys, new = new, old = old)
    } else stop("differing lengths")

    return(rtn)
}
compare_qrys <- function(x, new, old) {
    remove_spaces(old[[x]]) == remove_spaces(new[[x]]$query)
}
remove_spaces <- function(x) str_replace_all(x, " ", "")
 
# getargs <- list("validated" = list(table = "validated", mirna = "hsa-miR-199a-3p"),
#                 "mir2disease" = list(table = "mir2disease", mirna = "hsa-miR-199a-3p", 
#                                      target = "TP53", org = "rno", disease.drug = "cisplatin"),
#                 "disease.drug" = list(table = "disease.drug", mirna = "hsa-miR-199a-3p",
#                                       #target = "TP53",
#                                       org = "rno", disease.drug = "cisplatin"),
#                 "validated" = list(table = "predicted", mirna = "hsa-miR-199a-3p",
#                                    target = "TP53", org = "rno", predicted.cutoff.type = "p",
#                                    predicted.cutoff = 70, predicted.site = "all"))
# Test org parsing
arg_test_org <- 
    data_frame(table = "validated", 
               mirna = "hsa-miR-199a-3p", 
               org   = c("hsa", "human", "homo sapiens", "mmu",
                         "mouse", "mus musculus", "rno", "rat", 
                         "rattus norvegicus")) %>%
    test_arg_df

# Test disease.drug case-insensitive
arg_test_disease.drug <- 
    expand.grid(table = "disease.drug", 
                org = "rno", 
                disease.drug = c("bladder cancer", "cisplatin", "CiSplatin",
                                 "blaDDER Cancer"),
                stringsAsFactors = FALSE) %>% 
    test_arg_df

# Test multiple disease.drug, mirna, or target
arg_test_multiples <-
    list(
         data_frame(table = "all", 
                    disease.drug = list(c("bladder cancer", "cisplatin"))), 
         data_frame(table = "all",
                    mirna = list(c("MIMAT0000072", "hsa-miR-199a-3p")))
         ) %>%
    map(test_arg_df)
# checkreturneddata <-
#     map(arg_test_multiples$queries[[1]], ~ map(.x, search.multimir)) %>%
#     pmap(~ identical(.x, .y))


old <- get.multimirold(table = "disease.drug", 
                       #org = c("mmu", "hsa", "rno"),
                       disease.drug = c("cisplatin", "bladder cancer")



qry_disease     <- get.multimir(table = "mir2disease", mirna = "hsa-miR-199a-3p",
                                target = "TP53", org = "rno", 
                                disease.drug = "cisplatin")

list(data.frame(table = "predicted", mirna = "hsa-miR-199a-3p",
                target = "TP53", org = "rno", predicted.cutoff.type = "p",
                predicted.cutoff = 70, predicted.site = "all")) %>%
map(test_arg_df)

new_pred <- get.multimir(org = "mmu", 
                         target  = "Gnb1",
                         table   = "targetscan",
                         #summary = TRUE,
                         predicted.cutoff      = 35,
                         predicted.cutoff.type = "p",
                         predicted.site        = "all")
old_pred <- get.multimirold(org = "mmu", 
                            target  = "Gnb1",
                            table   = "targetscan",
                            #summary = TRUE,
                            predicted.cutoff      = 35,
                            predicted.cutoff.type = "p",
                            predicted.site        = "conserved")



