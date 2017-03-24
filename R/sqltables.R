#
# SQL query features
# 
# The queries have 4 key parts: SELECT, FROM, ON, and WHERE.
#       For each query characteristic:
#       .select is a single string of variables separated by commas
#       .from is a vector of strings, each of the format "table AS initial"
# The database has 3 mirna 'types': Validated, Predicted, and Associations, plus
#   3 additional filters...
#
# Table types
validated_tables <- function() c("mirecords", "mirtarbase", "tarbase")
predicted_tables <- function() c("diana_microt", "elmmo", "microcosm",
                                 "miranda", "mirdb", "pictar", "pita",
                                 "targetscan")
associations_tables <- function() c("mir2disease", "pharmaco_mir", "phenomir")
# Other key categories
tables_wo_target    <- function() c("mir2disease", "phenomir")
conserved_tables    <- function() c("miranda", "pita", "targetscan")


sql_validated <- function(.table) {

    this_type <- .table %in% validated_tables()
    as_mmsql(.select = if (this_type) c("i.experiment, i.support_type, i.pubmed_id")  else NULL,
             .from   = if (this_type) "%s AS i" else NULL)

}

sql_predicted <- function(.table) {

    this_type <- .table %in% predicted_tables()
    score_var <- get_score_var(.table)

    as_mmsql(.select = if (!this_type) NULL else sprintf("%s AS score", score_var),
             .from   = if (!this_type) NULL else sprintf("%s AS i", .table))

}

sql_associations <- function(.table) {
    
    this_type <- .table %in% associations_tables()

    # Build select list
    assoc   <- switch(.table, pharmaco_mir = "i.drug", "i.disease")
    pubmed  <- switch(.table,
                      mir2disease = "CONCAT_WS('. ', i.year, i.title)",
                      "i.pubmed_id")
    .select <- sprintf("%s AS disease_drug, %s AS paper_pubmedID", assoc, pubmed)
    # Build where var list
    wherevar <- switch(.table,
                       pharmaco_mir = "i.drug",
                       mir2disease = "i.disease", 
                       phenomir = "i.disease", "i.disease_class")
    .where  <- as_where(.vars = wherevar, .connect  = "OR", .operator = "IN")

    as_mmsql(.select     = if (!this_type) NULL else .select,
             .from       = if (!this_type) NULL else sprintf("%s AS i", .table),
             .where_list = as_where_list(.where))

}

sql_mirna <- function() {

    .where <- as_where(.vars     = c("m.mature_mirna_acc", "m.mature_mirna_id"),
                       .connect  = "OR",
                       .operator = "IN")
    as_mmsql(.select     = c("m.mature_mirna_acc, m.mature_mirna_id"),
             .from       = "mirna AS m",
             .on         = "m.mature_mirna_uid = i.mature_mirna_uid",
             .where_list = as_where_list(.where))

}

sql_target <- function(.table) {

    no_target <- .table %in% tables_wo_target()
    vars      <- c("t.target_symbol", "t.target_entrez", "t.target_ensembl")
    .select   <- paste(paste0("%1$s", vars), collapse = ", ")
    na_txt    <- "'NA' AS "
    .where    <- as_where(.vars     = if (no_target) NULL else vars,
                          .connect  = "OR",
                          .operator = "IN")

    as_mmsql(.select = if (no_target) sprintf(.select, na_txt) else sprintf(.select, ""),
             .from   = if (no_target) NULL else "target AS t",
             .on     = if (no_target) NULL else "i.target_uid = t.target_uid",
             .where_list = as_where_list(.where))

}

