################################################################################
# Script for testing during active development
################################################################################

# Load packages
library(devtools)
#library(roxygen2)
library(dplyr)
library(tidyr)
library(purrr)
#library(Rd2roxygen)


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
install.packages("test/multiMiRold_0.98.0.1.tar.gz")
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
qry_disease <- get.multimir(table = "disease.drug", mirna = "hsa-miR-199a-3p", #target = "TP53",
                            org = "rno", disease.drug = "cisplatin")


