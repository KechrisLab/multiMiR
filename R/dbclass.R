library(stringr)

# Create and assign class in one step
dbtable <- structure(list(), class = "")


# dbtable 
dbtable <- function(type = c("validated", "predicted", "disease", "drug"), 
                    columns = NULL) {

    type <- match.arg(type)
    has_target <- ifelse(type %in% 

    other_selects <- c("score_txt", "disease_txt", "pubmed_txt")
    
    .class <- sprintf("mm%s", type)

    score <- if (type == 'predicted') .score else NULL 

    # type dependent list items
#     list(score,
#          has_target,

    function(table_name, columns = NULL, has_target) {
        structure(c(table_name,
                    .target(table_name),
                    .score(table_name)

                       ),
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

mmquery <- function(table_type, 

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






# Defines columns to select for a dbtable object
# 
# @param x vector of additional columns to select
# @keywords internal
.col_select <- function(x = NULL) {
    c(c("m.mature_mirna_acc", "m.mature_mirna_id", 
        "%s t.target_symbol", "%s t.target_entrez",
        "%s t.target_ensembl"), x) 
}

expand_select <- function(.col_select) {
    paste(.col_select, collapse = ", ") 
}


.score <- function(table_name) {
    switch(table_name,
           diana_microt = "i.miTG_score",
           elmmo        = "i.p",
           microcosm    = "i.score",
           mirdb        = "i.score",
           pictar       = "i.score",
           miranda      = "i.mirsvr_score",
           pita         = "i.ddG",
           targetscan  =  "i.context_plus_score")
}

# Feature of all
.target <- function(table_name) {

    no_target <- c("mir2disease", "phenomir")
    if (table_name %in% no_target) {
        .select = "'NA' AS"
        .uid = ""
        where_target = FALSE
    } else {
        .select = ""
        .uid = "AND i.target_uid=t.target_uid"
        where_target = TRUE
    }

    return(list(.select = .select, .uid = .uid, where_target))
}

.diseasedrug <- function(table_name) {
    .select <- switch(table_name,
                      pharmaco_mir = "i.drug",
                      "i.disease")
    .pubmed <- switch(table_name,
                      mir2disease = "CONCAT_WS('. ', i.year, i.title)",
                      "i.pubmed_id")
    .query  <- switch(table_name, 
                      mir2disease  = structure(c("i.disease IN"),class = "where_or"),
                      pharmaco_mir = structure(c("i.drug IN"), class = "where_or"),
                      phenomir     = structure(c("i.disease IN", "i.disease_class IN"),
                                               class = "sqlwhere", "in", "or"))
    
    list(.select = .select, .pubmed = .pubmed, .query = .query)
}





# wrap each WHERE arg in quotes, then all of them in parentheses
where_value_in <- function(values) {
    paste0("(", paste( paste0("'", values, "'"), collapse = ", "), ")")
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




