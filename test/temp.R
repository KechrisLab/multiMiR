install.packages("RMySQL")
library(DBI)
library(RMySQL)
library(dplyr)

mm_url <- "http://multimir.ucdenver.edu/multiMiR_dbTables.txt"

RMySQL::MySQL()
con <- dbConnect

src_mysql(dbname = "multimir", host = "http://multimir.ucdenver.edu/cgi-bin/multimir_univ.pl") 

tbl(

