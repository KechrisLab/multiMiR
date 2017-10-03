#' Functions defining the category each table belongs to.  
#' 
#' One of three types: predicted, validated, or diseasedrug.
#' Additionally two functions define characteristics of tables: those without a
#' target column \code{tables_wo_target} and those with conserved target sites
#' \code{conserved_tables}.
#'
#' @param .table a table name
#' @return Returns dataset that names that belong to the category of the
#' function name (e.g. \code{validated_tables()} returns tables with validated
#' miRNA-target interactions).  \code{reverse_table_lookup()} does the opposite;
#' it returns the category a given \code{.table} belongs to.
#'
#' @examples
#' all_tables()
#' validated_tables()
#' predicted_tables()
#' diseasedrug_tables()
#' predicted_tables() %in% all_tables() # TRUE
#' table_types()
#'
#' @aliases all_tables, validated_tables, predicted_tables, diseasedrug_tables,
#' tables_wo_target, reverse_table_lookup
#' @keywords tables
#' @export
all_tables <- function() {
    c(validated_tables(), predicted_tables(), diseasedrug_tables())
}

#' @rdname all_tables
#' @export
validated_tables <- function() c("mirecords", "mirtarbase", "tarbase")

#' @rdname all_tables
#' @export
predicted_tables <- function() c("diana_microt", "elmmo", "microcosm",
                                 "miranda", "mirdb", "pictar", "pita",
                                 "targetscan")

#' @rdname all_tables
#' @export
diseasedrug_tables <- function() c("mir2disease", "pharmaco_mir", "phenomir")

#' @rdname all_tables
#' @export
tables_wo_target    <- function() c("mir2disease", "phenomir")

#' @rdname all_tables
#' @export
conserved_tables    <- function() c("miranda", "pita", "targetscan")

#' @rdname all_tables
#' @export
reverse_table_lookup <- function(.table) {
    if (.table %in% validated_tables()) {
        "validated"
    } else if (.table %in% predicted_tables()) {
        "predicted"
    } else if (.table %in% diseasedrug_tables()) {
        "disease.drug"
    }
}

#' @rdname all_tables
#' @export
table_types <- function() c("validated", "predicted", "disease.drug")

#' @keywords internal
get_score_var <- function(.table) {
    switch(.table,
           diana_microt = "i.miTG_score",
           elmmo        = "i.p",
           microcosm    = "i.score",
           mirdb        = "i.score",
           pictar       = "i.score",
           miranda      = "i.mirsvr_score",
           pita         = "i.ddG",
           targetscan  =  "i.context_plus_score",
           NULL)
}

