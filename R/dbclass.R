
# Create and assign class in one step
dbtable <- structure(list(), class = "")


#' dbtable 
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



predicted_table <- function(table_name, columns = NULL, has_target) {
    structure(list(table_name,
                   .col_select(columns),
                   has_target,
                   score(table_name)
                   ), 
              class = "mmpredicted")
}

print

def_dbtable <- function() {

    UseMethod("def_dbtable")

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
           "i.miTG_score"         = "diana_microt",
           "i.p "                 = "elmmo",
           "i.score"              = "microcosm",
           "i.score"              = "mirdb",
           "i.score"              = "pictar",
           "i.mirsvr_score"       = "miranda",
           "i.ddG"                = "pita",
           "i.context_plus_score" = "targetscan") 
}

.target <- function(table_name) {

    no_target <- c("mir2disease", "phenomir")
    if (table_name %in% no_target) {
        .select = "'NA' AS"
        .uid = ""
    } else {
        .select = ""
        .uid = "AND i.target_uid=t.target_uid"
    }

    return(c(.select = .select, .uid    = .uid))
}

