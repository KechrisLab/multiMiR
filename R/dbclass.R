{ 
library(stringr)

# Create and assign class in one step
dbtable <- structure(list(), class = "")


# dbtable 
dbtable <- function(type = c("validated", "predicted", "disease", "drug"), 
                    columns = NULL) {

    type <- match.arg(type)
#     has_target <- ifelse(type %in% 

    other_selects <- c("score_txt", "disease_txt", "pubmed_txt")
    
    .class <- sprintf("mm%s", type)

    score <- if (type == 'predicted') .score else NULL 

    # type dependent list items
#     list(score,
#          has_target,

    function(table_name, columns = NULL, has_target) {
        structure(c(table_name,
                    .target(table_name),
                    .score(table_name)),
                  class = .class)
    }
}

# org = "hsa", 
# mirna = NULL, 
# target = NULL,
# disease.drug = NULL, 
# table = "validated",
# predicted.cutoff = NULL,
# predicted.cutoff.type = "p",
# predicted.site = "conserved", 

predtable    <- structure(list(argslist = NULL, target = TRUE, score = TRUE, disease_drug = FALSE), class = "mmpredicted")
validtable   <- structure(list(argslist = NULL, target = TRUE, score = FALSE, disease_drug = FALSE), class = "mmpredicted")
diseasetable <- structure(list(argslist = NULL, target = TRUE, score = FALSE, disease_drug = TRUE), class = "mmpredicted")


build_query <- function(mmtable_obj) {


}


predicted_table <- function(table_name, columns = NULL, has_target) {
    structure(list(table_name,
                   .col_select(columns),
                   has_target,
                   score(table_name)
                   ), 
              class = "mmpredicted")
}






############################################################
# Structure by table features.  
#   For each new table, the question is now: does it have each of these features
#   and if so what are the select vars, from table and indicator, 'on' ids, and
#   where variables
############################################################

sql_feature <- function() {
    .select
    .from
    .on
    .where
}


collapse_vec <- function(x, sep = ", ") {

}

# wrap each WHERE arg in quotes, then all of them in parentheses
where_value_in <- function(values) {
    wrap_parens(paste( paste0("'", values, "'"), collapse = ", "))
}

# Combine and collapse WHERE ... IN (...) OR statements....
# .vars -> c("i.disease", "i.disease_class")
# .values -> c("hsa", "mmu")
# wherein_vars(.vars, wherein_value(.values)) == "i.disease IN ('hsa', 'mmu') # OR i.disease_class IN ('hsa', 'mmu')"
where_vars <- function(vars, value, type = "IN") {
    stopifnot(str_detect(type, "IN|[><=]"))
    paste(paste(vars, value, sep = pad(type)), collapse = " OR ")
}

sql_where <- function(vars, values, type = "IN") {

    type 
    switch(type,
           IN   = where_vars(vars, where_value_in(values)),
           `>=` = where_vars(vars, values, type = ">="),
           `=`  = where_vars(vars, values, type = "="))

}

pad <- function(x) paste0(" ", x, " ")

} 



this_table <- "mir2disease" %>% query_features
select_list <- this_table %>% transpose %>% .$.select
from_list <- this_table %>% transpose %>% .$.from


############################################################
# SQL query features
# 
# The queries have 4 key parts: SELECT, FROM, ON, and WHERE.
#       For each query characteristic:
#       .select is a single string of variables separated by commas
#       .from is a vector of strings, each of the format "table AS initial"
# The database has 3 mirna 'types': Validated, Predicted, and Associations, plus
#   3 additional filters...
############################################################
require(purrr)

validated <- function(table_name) {
    if (table_name %in% c("mirecords", "mirtarbase", "tarbase")) {
       .select <- c("i.experiment, i.support_type, i.pubmed_id") 
    } else {
        .select <- NULL
    }
    .from   <- NULL
    .on     <- NULL
    .where  <- NULL

    return(list(.select = .select, .from = .from, .on = .on, .where = .where))
}

predicted <- function(table_name) {
    score_var <- switch(table_name,
                        diana_microt = "i.miTG_score",
                        elmmo        = "i.p",
                        microcosm    = "i.score",
                        mirdb        = "i.score",
                        pictar       = "i.score",
                        miranda      = "i.mirsvr_score",
                        pita         = "i.ddG",
                        targetscan  =  "i.context_plus_score",
                        NULL)
    .select <- if (is.null(score_var)) NULL else sprintf("%s AS score", score_var)
    .from   <- NULL
    .on     <- NULL
    .where  <- TRUE

    return(list(.select = .select, .from = .from, .on = .on, .where = .where))
}

associations <- function(table_name) {
    
    assoc   <- switch(table_name,
                      pharmaco_mir = "i.drug",
                      "i.disease")
    pubmed  <- switch(table_name,
                      mir2disease = "CONCAT_WS('. ', i.year, i.title)",
                      "i.pubmed_id")
    .select <- sprintf("%s AS disease_drug, %s AS paper_pubmedID", assoc, pubmed)
    .from   <- NULL
    .on     <- NULL
    .where  <- TRUE

    return(list(.select = .select, .from = .from, .on = .on, .where = .where))
}

mirna <- function(table_name) {
    .select <- c("m.mature_mirna_acc, m.mature_mirna_id")
    .from   <- "mirna AS m"
    .on     <- "m.mature_mirna_uid = i.mature_mirna_uid"
    .where  <- c("m.mature_mirna_acc", "m.mature_mirna_id")
    return(list(.select = .select, .from = .from, .on = .on, .where = .where))
}

org <- function(table_name, target_where) {
    .select <- NULL
    .from   <- NULL
    .on     <- NULL
    .where  <- ifelse(is.null(target_where), c("m.org", "t.org"), "m.org")

    return(list(.select = .select, .from = .from, .on = .on, .where = .where))
}

target <- function(table_name) {
    no_target <- c("mir2disease", "phenomir")

    .select <- c("%st.target_symbol, %st.target_entrez, %st.target_ensembl")
    na_txt <- "'NA' AS "

    if (table_name %in% no_target) {
        .select <- sprintf(.select, na_txt, na_txt, na_txt)
        .from   <- NULL
        .on     <- NULL
        .where  <- NULL
    } else {
        .select <- sprintf(.select, "", "")
        .from   <- "target AS t"
        .on     <- "i.target_uid = t.target_uid"
        .where  <- c("t.target_symbol", "t.target_entrez", "t.target_ensembl")
    }


    return(list(.select = .select, .from = .from, .on = .on, .where = .where))
}



query_features <- function(table_name){ 
        validated    <- validated(table_name)
        predicted    <- predicted(table_name)
        associations <- associations(table_name)
        mirna        <- mirna(table_name)
        target       <- target(table_name)
        org          <- org(table_name, target$.where)

        rtn <- list(this_table   = table_name, 
                    mirna        = mirna,
                    target       = target,
                    validated    = validated,
                    predicted    = predicted,
                    associations = associations,
                    org          = org)
    }


expand_select <- function(select_list) {
    paste("SELECT", paste(unlist(select_list), collapse = ", "))
}

expand_from <- function(table_name, from_list) {
# NOTE: issue here is where to get the table_name from... Do we add it to the
    # query_features list as a separate item (as i did for now), or should it be
    # defined in validated(), predicted(), and associations()???
    # Kind of depends on what happens with the other queries: query_predicted()
    # etc. and the table lookup
}

expand_in <- function(in_list) {

}

