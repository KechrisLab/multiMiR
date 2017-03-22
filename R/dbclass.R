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



this_table  <- 
    "diana_microt" %>% 
    #"mir2disease" %>% 
    query_features
select_list <- this_table %>% transpose %>% .$.select
from_list   <- this_table %>% transpose %>% .$.from
on_list     <- this_table %>% transpose %>% .$.on

paste(expand_select(select_list),
      expand_from(from_list),
      expand_on(on_list))

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

validated <- function(.table) {

    these_tables <- c("mirecords", "mirtarbase", "tarbase")
    this_type <- .table %in% these_tables

    .select <- if (this_type) c("i.experiment, i.support_type, i.pubmed_id")  else NULL
    .from   <- if (this_type) "%s AS i" else NULL
    .on     <- NULL
    .where  <- NULL

    return(list(.select = .select, .from = .from, .on = .on, .where = .where))
}

predicted <- function(.table) {

    these_tables <- c("diana_microt", "elmmo", "microcosm", "miranda",
                      "mirdb", "pictar", "pita", "targetscan")
    this_type <- .table %in% these_tables

    score_var <- switch(.table,
                        diana_microt = "i.miTG_score",
                        elmmo        = "i.p",
                        microcosm    = "i.score",
                        mirdb        = "i.score",
                        pictar       = "i.score",
                        miranda      = "i.mirsvr_score",
                        pita         = "i.ddG",
                        targetscan  =  "i.context_plus_score")

    .select <- if (this_type) sprintf("%s AS score", score_var) else NULL
    .from   <- if (this_type) sprintf("%s AS i", .table) else NULL
    .on     <- NULL
    .where  <- TRUE

    return(list(.select = .select, .from = .from, .on = .on, .where = .where))
}

associations <- function(.table) {
    
    these_tables <- c("mir2disease", "pharmaco_mir", "phenomir")
    this_type <- .table %in% these_tables

    assoc   <- switch(.table, pharmaco_mir = "i.drug", "i.disease")
    pubmed  <- switch(.table,
                      mir2disease = "CONCAT_WS('. ', i.year, i.title)",
                      "i.pubmed_id")
    .select <- sprintf("%s AS disease_drug, %s AS paper_pubmedID", assoc, pubmed)

    .select <- if (!this_type) NULL else .select
    .from   <- if (!this_type) NULL else sprintf("%s AS i", .table)
    .on     <- NULL
    .where  <- TRUE

    return(list(.select = .select, .from = .from, .on = .on, .where = .where))
}

mirna <- function() {
    .select <- c("m.mature_mirna_acc, m.mature_mirna_id")
    .from   <- "mirna AS m"
    .on     <- "m.mature_mirna_uid = i.mature_mirna_uid"
    .where  <- c("m.mature_mirna_acc", "m.mature_mirna_id")
    return(list(.select = .select, .from = .from, .on = .on, .where = .where))
}

org <- function(target_where) {
    .select <- NULL
    .from   <- NULL
    .on     <- NULL
    .where  <- ifelse(is.null(target_where), c("m.org", "t.org"), "m.org")

    return(list(.select = .select, .from = .from, .on = .on, .where = .where))
}

target <- function(.table) {

    tables_wo_target <- c("mir2disease", "phenomir")
    no_target        <- .table %in% tables_wo_target

    vars    <- c("t.target_symbol", "t.target_entrez", "t.target_ensembl")
    .select <- paste(paste0("%1$s", vars), collapse = ", ")
    na_txt  <- "'NA' AS "

    .select <- if (no_target) sprintf(.select, na_txt) else sprintf(.select, "")
    .from   <- if (no_target) NULL else "target AS t"
    .on     <- if (no_target) NULL else "i.target_uid = t.target_uid"
    .where  <- if (no_target) NULL else vars

    return(list(.select = .select, .from = .from, .on = .on, .where = .where))
}



query_features <- function(.table){ 
        validated    <- validated(.table)
        predicted    <- predicted(.table)
        associations <- associations(.table)
        mirna        <- mirna()
        target       <- target(.table)
        org          <- org(target$.where)

        rtn <- list(#this_table   = .table, 
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

expand_from <- function(from_list) {

    from_list <- merge_order(from_list)
    paste("FROM", paste(unlist(from_list), collapse = ", "))

}

expand_on <- function(on_list) {

    on_list <- merge_order(on_list)
    paste0("ON (", paste(unlist(on_list), collapse = " AND "), ")")

}

merge_order <- function(.list) {
    c(.list["mirna"],
      .list[c("validated", "predicted", "associations")],
      .list["target"])
}

