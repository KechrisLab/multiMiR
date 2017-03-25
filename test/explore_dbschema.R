



########################################
# Test/explore new s3 approach 

this_table  <- 
    #"diana_microt" %>% 
    #"mir2disease" %>% 
    #"targetscan" %>% 
sql_features(.table = "targetscan", org, predicted.site, 
             predicted.cutoff.type = "p", 
             predicted.cutoff = 35)

select_list <- this_table %>% transpose %>% .$.select
from_list   <- this_table %>% transpose %>% .$.from
on_list     <- this_table %>% transpose %>% .$.on
where_list  <- this_table %>% transpose %>% .$.where

paste(expand_select(select_list),
      expand_from(from_list),
      expand_on(on_list))

combine_wheres(where_list)



example3 <- get.multimir(org = "mmu", 
                         target  = "Gnb1",
                         table   = "miranda",
                         summary = TRUE,
                         predicted.cutoff      = 35,
                         predicted.cutoff.type = "p",
                         predicted.site        = "all")
org <- "mmu" 
target  <- "Gnb1"
table   <- "miranda"
summary <- TRUE
predicted.cutoff      <- 35
predicted.cutoff.type <- "p"
predicted.site        <- "all"

build_mmquery(.table = table, org = "mmu", mirna = NULL, target = target, 
              predicted.cutoff      = predicted.cutoff, 
              predicted.cutoff.type = predicted.cutoff.type,
              predicted.site        = predicted.site)

build_mmquery(.table = "targetscan")






















################################################################################











































# VALIDATED: 
#   base + (mirna and/or target) + org
# DISEASE: Basic SQL structure is:
#   mir2disease:  base + disease.drug + org + mirna
#   pharmaco_mir: base + target + disease.drug + org + mirna
#   phenomir:     base + disease.drug + org + mirna
# PREDICTED: Basic SQL structure:
#   (choose score var) + base + (mirna and/or target) + org + 
#       conserved(for 3 of 8) + cutoff

library(multiMiR)
library(dplyr)
library(tidyr)
library(pryr)
library(purrr)
library(stringr)
library(forcats)

# target genes table (target)
"SELECT * FROM target LIMIT 20" %>% search.multimir
"SELECT * FROM target WHERE org = 'rno' LIMIT 20" %>% search.multimir
"SELECT * FROM target WHERE org = 'mmu' LIMIT 20" %>% search.multimir

# microRNA table (mirna)
"SELECT * FROM mirna LIMIT 20" %>% search.multimir #hsa
"SELECT * FROM mirna WHERE org = 'rno' LIMIT 20" %>% search.multimir
"SELECT * FROM mirna WHERE org = 'mmu' LIMIT 20" %>% search.multimir


validated    = c("mirecords", "mirtarbase", "tarbase"),
predicted    = c("diana_microt", "elmmo", "microcosm",
                 "miranda", "mirdb", "pictar", "pita",
                 "targetscan"),
disease.drug = c("mir2disease", "pharmaco_mir", "phenomir"),


# Validated tables
"SELECT * FROM mirecords LIMIT 20" %>% search.multimir %>% tbl_df


schema <- 
    readLines(multiMiR:::full_url("multimir.schema")) 

dbguide <- 
    data_frame(text = schema) %>% 
        filter(!str_detect(text, "^--")) %>% 
        filter(text != "") %>% 
        mutate(text = str_replace_all(text, "[\t]+", "")) %>%
        separate(text, c("sql", "note"), "--", fill = "right", remove = TRUE) %>%
        mutate(table = sql %>% str_extract("(?<=\\`)(.*?)(?=\\`)")) %>%
        fill(table) 

dbguide %>% filter(str_detect(note, "mature"))

dbguide %>% 
    mutate(DROP = str_detect(sql, "DROP")) %>%
    mutate(CREATE = str_detect(sql, "CREATE")) %>%
    mutate(AUTO_INCREMENT = str_detect(sql, "AUTO_INCREMENT")) %>%
    mutate(KEY = str_detect(sql, "KEY") & !str_detect(sql, "PRIMARY"))

# validated's wheres: mirnatarget, org
# predicted's wheres: mirnatarget, org, conserved, cutoff
# diseased's  wheres: mirnatarget, org, diseasedrug

paste("SELECT m.mature_mirna_acc, m.mature_mirna_id,",
      "t.target_symbol, t.target_entrez, t.target_ensembl,",
      "i.experiment, i.support_type, i.pubmed_id",
      "FROM", mirna.table, "AS m INNER JOIN", table, 
      "AS i INNER JOIN", target.table, 
      "AS t ON (m.mature_mirna_uid=i.mature_mirna_uid AND",
      "i.target_uid=t.target_uid) WHERE")


# predicted
paste("SELECT m.mature_mirna_acc, m.mature_mirna_id,",
      "t.target_symbol, t.target_entrez, t.target_ensembl,",
      score_vars, "AS score", 
      "FROM", mirna.table, "AS m INNER JOIN", table, 
      "AS i INNER JOIN", target.table, 
      "AS t ON (m.mature_mirna_uid=i.mature_mirna_uid",
      "AND i.target_uid=t.target_uid) WHERE")
# mir2disease 
paste("SELECT m.mature_mirna_acc, m.mature_mirna_id, 'NA' AS",
      "target_symbol, 'NA' AS target_entrez, 'NA' AS",
      "target_ensembl, i.disease AS disease_drug,",
      "CONCAT_WS('. ', i.year, i.title) AS paper_pubmedID",
      "FROM", mirna.table, "AS m INNER JOIN", table, 
      "AS i ON (m.mature_mirna_uid=i.mature_mirna_uid) WHERE")
# pharmaco_mir 
paste("SELECT m.mature_mirna_acc, m.mature_mirna_id,",
      "t.target_symbol, t.target_entrez, t.target_ensembl,",
      "i.drug AS disease_drug, i.pubmed_id AS paper_pubmedID",
      "FROM", mirna.table, "AS m INNER JOIN", table, 
      "AS i INNER JOIN", target.table, 
      "AS t ON (m.mature_mirna_uid=i.mature_mirna_uid AND",
      "i.target_uid=t.target_uid) WHERE")
# phenomir     
paste("SELECT m.mature_mirna_acc, m.mature_mirna_id, 'NA' AS",
      "target_symbol, 'NA' AS target_entrez, 'NA' AS",
      "target_ensembl, i.disease AS disease_drug, i.pubmed_id AS",
      "paper_pubmedID FROM", mirna.table, "AS m INNER JOIN", table, 
      "AS i ON (m.mature_mirna_uid=i.mature_mirna_uid) WHERE")




# predicted/validated
# SELECT m.mature_mirna_acc, m.mature_mirna_id, t.target_symbol, t.target_entrez,
#       t.target_ensembl, 
# (predicted only: score_vars; 'AS score' for some tables)
# FROM mirna.table AS m INNER JOIN table AS i INNER JOIN target.table AS t ON
# (m.mature_mirna_uid=i.mature_mirna_uid", "AND i.target_uid=t.target_uid)

# then compared with disease.drug
# SELECT m.mature_mirna_acc, m.mature_mirna_id, 
#   t.target_symbol, t.target_entrez, t.target_ensembl, # =NA for mir2disease and phenomir
# (disease only: 
#      mir2disease: 1) i.disease AS disease_drug 
#                   2) CONCAT_WS('. ', i.year, i.title) AS paper_pubmedID
#      pharmaco_mir: 1) i.drug AS disease_drug 
#                    2) i.pubmed_id AS paper_pubmedID
#      phenomir:    1) i.disease AS disease_drug
#                   2) i.pubmed_id AS paper_pubmedID
# FROM mirna.table AS m INNER JOIN table AS i INNER JOIN target.table AS t ON
# (m.mature_mirna_uid=i.mature_mirna_uid", "AND i.target_uid=t.target_uid) [target only on pharmaco_mir]


SELECT <- c("m.mature_mirna_acc", "m.mature_mirna_id", "t.target_symbol", "t.target_entrez")

gen_query_metadata <- function(table_name, type, has_target, score_var,
                               as_score, disease_drug_var, pubmed_var) {

    if (type == 'predicted' & (is.null(score_var) | is.null(as_score))) { 
        stop('score column name and as_score=T/F required for predicted tables')
    }
    if (type == "disease_drug" & (is.null(disease_drug) | is.null(paper_pubmedID))) {
        stop('disease_drug and pubmed_var required for disease_drug tables')
    }

    data_frame(table_name       = table_name,
               type             = type,
               has_target       = has_target,
               score_var        = score_var,
               as_score         = as_score,
               disease_drug_var = disease_drug_var,
               pubmed_var       = pubmed_var
               )

}

# target: mir2disease, phenomir == FALSE, else TRUE
# predicted: score_vars == list score var for each, 
#   predicted 'AS score': microcosm, mirdb, pictar == FALSE, else TRUE
# disease_drug, disease_drug == list disease_drug var for each
#               paper_pubmedID == list pubmedid var for each

# na's: if target=FALSE then 'NA' AS t.vars
# ON: if target= TRUE then "AND i.target_uid=t.target_uid"

# create_basemmdf <- function(type, table_name) {
#     tibble::data_frame(type = type,
#                        table_name = table_name) 
# }
# list(
# 
# list(
#      create_basemmdf(type       = "validated",
#                      table_name = c("mirecords", "mirtarbase",
#                                     "tarbase")),
#      create_basemmdf(type       = "predicted",
#                      table_name = c("diana_microt", "elmmo",
#                                     "microcosm", "miranda", "mirdb",
#                                     "pictar", "pita", "targetscan")),
#      create_basemmdf(type       = "predicted",
#                      table_name = c("mir2disease", "pharmaco_mir",
#                                     "phenomir"))
#      )

# mmdb_metadata %>%
query_metadata <- 
    data_frame(type = c(rep("validated", 3), rep("predicted", 8), rep("disease.drug", 3)),
               table_name   = c("mirecords", "mirtarbase", "tarbase",
                                "diana_microt", "elmmo", "microcosm",
                                "miranda", "mirdb", "pictar", "pita",
                                "targetscan",
                                "mir2disease", "pharmaco_mir", "phenomir")
               ) %>%
    mutate(has_target = if_else(table_name %in% c("mir2disease", "phenomir"), FALSE, TRUE),
           score_var  = table_name %>%
               fct_recode("i.miTG_score"         = "diana_microt",
                          "i.p "                 = "elmmo",
                          "i.score"              = "microcosm",
                          "i.score"              = "mirdb",
                          "i.score"              = "pictar",
                          "i.mirsvr_score"       = "miranda",
                          "i.ddG"                = "pita",
                          "i.context_plus_score" = "targetscan") %>% as.character,
           score_var = score_var %>% ifelse(. == table_name, NA, .),
           as_score = ifelse(type == "predicted", TRUE, NA) %>%
                      if_else(table_name %in% c("microcosm", "mirdb", "pictar") & 
                              type == "predicted", FALSE, .)) 
# Add disease_drug info
query_metadata <- 
    query_metadata %>%
    mutate(disease_drug = table_name %>%
               fct_recode("i.disease" = "mir2disease", 
                          "i.disease" = "phenomir",
                          "i.drug"    = "pharmaco_mir") %>% as.character,
           disease_drug = disease_drug %>% ifelse(. == table_name, NA, .)) %>%
    mutate(pubmed_var = table_name %>%
               fct_recode("CONCAT_WS('. ', i.year, i.title)" = "mir2disease", 
                          "i.pubmed_id" = "phenomir",
                          "i.pubmed_id" = "pharmaco_mir") %>% as.character,
           pubmed_var = pubmed_var %>% ifelse(. == table_name, NA, .)) 
# Add WHERE info
query_metadata <- 
    query_metadata %>%
    mutate(where_mirna = TRUE,
           where_target = has_target,
           where_org = TRUE,
           where_conserved = ifelse(type == "predicted", TRUE, FALSE),
           where_cutoff    = where_conserved, 
           where_disease   = ifelse(type == "disease.drug", TRUE, FALSE))



#      mir2disease: 1) i.disease AS disease_drug 
#                   2) CONCAT_WS('. ', i.year, i.title) AS paper_pubmedID
#      pharmaco_mir: 1) i.drug AS disease_drug 
#                    2) i.pubmed_id AS paper_pubmedID
#      phenomir:    1) i.disease AS disease_drug
#                   2) i.pubmed_id AS paper_pubmedID
# validated's wheres: mirnatarget, org
# predicted's wheres: mirnatarget, org, conserved, cutoff
# diseased's  wheres: mirnatarget, org, diseasedrug

# predicted/validated
# SELECT m.mature_mirna_acc, m.mature_mirna_id, t.target_symbol, t.target_entrez,
#       t.target_ensembl, 
# (predicted only: score_vars; 'AS score' for some tables)
# FROM mirna.table AS m INNER JOIN table AS i INNER JOIN target.table AS t ON
# (m.mature_mirna_uid=i.mature_mirna_uid", "AND i.target_uid=t.target_uid)

# then compared with disease.drug
# SELECT m.mature_mirna_acc, m.mature_mirna_id, 
#   t.target_symbol, t.target_entrez, t.target_ensembl, # =NA for mir2disease and phenomir
# (disease only: 
#      mir2disease: 1) i.disease AS disease_drug 
#                   2) CONCAT_WS('. ', i.year, i.title) AS paper_pubmedID
#      pharmaco_mir: 1) i.drug AS disease_drug 
#                    2) i.pubmed_id AS paper_pubmedID
#      phenomir:    1) i.disease AS disease_drug
#                   2) i.pubmed_id AS paper_pubmedID
# FROM mirna.table AS m INNER JOIN table AS i INNER JOIN target.table AS t ON
# (m.mature_mirna_uid=i.mature_mirna_uid", "AND i.target_uid=t.target_uid) [target only on pharmaco_mir]

gen_select <- function() {
    # score_txt, disease_txt, pubmed_txt need ', ' prefix if present
    # target_uid needs 'AND ' prefix if present
    sprintf("SELECT m.mature_mirna_acc, m.mature_mirna_id, 
                    %s t.target_symbol, %s t.target_entrez, %s t.target_ensembl 
                    %s %s %s
            FROM mirna.table AS m INNER JOIN table AS i INNER JOIN target.table AS t 
            ON (m.mature_mirna_uid=i.mature_mirna_uid %s) 
            WHERE ",
            target_na, target_na, target_na, score_txt, disease_txt, pubmed_txt,
            target_uid)
}


.target <- function(target = TRUE) {

    if (target) {
        .select = ""
        .uid = "AND i.target_uid=t.target_uid"
    } else {
        .select = "'NA' AS"
        .uid = ""
    }

    return(list(.select = .select,
                .uid    = .uid))
}

.score <- function(table_name) {
    table_name %>%
        switch("i.miTG_score"         = "diana_microt",
               "i.p "                 = "elmmo",
               "i.score"              = "microcosm",
               "i.score"              = "mirdb",
               "i.score"              = "pictar",
               "i.mirsvr_score"       = "miranda",
               "i.ddG"                = "pita",
               "i.context_plus_score" = "targetscan") %>% as.character
}


# cryptic score error: miranda, pita, targetscan
 


my_mirna <- 
    get.multimir(org = "mmu", 
                 target  = "Gnb1",
                 table   = "predicted",
                 summary = TRUE,
                 predicted.cutoff      = 35,
                 predicted.cutoff.type = "p",
                 predicted.site        = "all") %>% str
lapply( search.multimir )


get.multimir(org = "mmu", 
             target  = "Gnb1",
             table   = "diana_microt",
             summary = TRUE,
             predicted.cutoff      = 35,
             predicted.cutoff.type = "p",
             predicted.site        = "all") %>% .[[1]] %>% search.multimir(query = .)
