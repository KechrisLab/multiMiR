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

table <- "miranda"
mirna <- "mirna"
target <- "target"
org <- "rno"
mirna.table  <- "mirna.table"
target.table <- "target.table"
predicted.cutoff <- NULL
predicted.cutoff.type <- "n"
predicted.site <- "conserved"

########################################
# Convert to roxygen2 documentation
########################################

# formatR::usage(Rd2roxygen)
# setwd("~/Projects/KechrisLab/")
# Rd2roxygen("multiMiR")

