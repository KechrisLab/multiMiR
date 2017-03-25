#' Generate mmsql objects for each of the three types of tables, as well as the
#' mirna and target tables.
#'
#' The three types of tables are predicted, validated, and associations
#' (disease/drug). 
#' 
#' @aliases sql, tables, table_types
#' @keywords tables types, predicted, validated, associations, disease, drug
#' 
#' @keywords internal
sql_validated <- function(.table) {

    this_type <- .table %in% validated_tables()
    as_mmsql(.select = if (!this_type) NULL else c("i.experiment, i.support_type, i.pubmed_id"),
             .from   = if (!this_type) NULL else sprintf("%s AS i", .table))

}

#' @rdname sql_validated
#' @keywords internal
sql_predicted <- function(.table, org, predicted.site, predicted.cutoff.type,
                          predicted.cutoff) {

    this_type <- .table %in% predicted_tables()

    if (!this_type) {
        .where_list <- NULL 
        .orderby <- NULL 
    } else {
        score_var    <- get_score_var(.table)
        cutoff_name  <- create_cutoff_name(.table, org, predicted.site)
        score_cutoff <- cutoff_to_score(.table, cutoff_name, predicted.cutoff.type,
                                        predicted.cutoff)
        conserved    <- where_conserved(.table, org, predicted.site)
        cutoff       <- where_cutoff(.table, score_var, score_cutoff)
        .orderby     <- as_orderby(.vars = score_var, .order = "DESC")
        .where_list  <- as_where_list(conserved = conserved,
                                      cutoff    = cutoff)
    }
    as_mmsql(.select = if (!this_type) NULL else sprintf("%s AS score", score_var),
             .from   = if (!this_type) NULL else sprintf("%s AS i", .table),
             .where_list  = .where_list,
             .orderby     = .orderby)

}

#' @rdname sql_validated
#' @keywords internal
sql_associations <- function(.table, disease.drug) {
    
    this_type <- .table %in% associations_tables()

    # Build select list
    assoc   <- switch(.table, pharmaco_mir = "i.drug", "i.disease")
    pubmed  <- switch(.table,
                      mir2disease = "CONCAT_WS('. ', i.year, i.title)",
                      "i.pubmed_id")
    .select <- sprintf("%s AS disease_drug, %s AS paper_pubmedID", assoc, pubmed)
    .where  <- where_associations(.table, disease.drug)

    as_mmsql(.select     = if (!this_type) NULL else .select,
             .from       = if (!this_type) NULL else sprintf("%s AS i", .table),
             .where_list = if (!this_type) NULL else as_where_list(.where))

}

#' @rdname sql_validated
#' @keywords internal
sql_mirna <- function(mirna) {

    if (is.null(mirna)) {
        .where_list <- NULL
    } else {
        .where <- as_where(.vars     = c("m.mature_mirna_acc", "m.mature_mirna_id"),
                           .connect  = "OR",
                           .operator = "IN",
                           .value    = mirna)
        .where_list <- as_where_list(.where)
    }
    as_mmsql(.select     = c("m.mature_mirna_acc, m.mature_mirna_id"),
             .from       = "mirna AS m",
             .on         = "m.mature_mirna_uid = i.mature_mirna_uid",
             .where_list = .where_list)

}

#' @rdname sql_validated
#' @keywords internal
sql_target <- function(.table, target) {

    no_target <- .table %in% tables_wo_target()
    prefix    <- if (no_target) "" else "t."
    na_txt    <- if (no_target) "'NA' AS " else ""
    vars      <- purrr::map_chr(c("%1$starget_symbol", "%1$starget_entrez",
                                  "%1$starget_ensembl"), sprintf, prefix)
    .select   <- paste(paste0("%1$s", vars), collapse = ", ")
    if (is.null(target)) {
        .where_list <- NULL
    } else {
        .where <- as_where(.vars     = if (no_target) NULL else vars,
                           .connect  = "OR",
                           .operator = "IN",
                           .value    = target)
        .where_list <- as_where_list(.where)
    }

    as_mmsql(.select = sprintf(.select, na_txt), #else sprintf(.select, "", ),
             .from   = if (no_target) NULL else "target AS t",
             .on     = if (no_target) NULL else "i.target_uid = t.target_uid",
             .where_list = .where_list)

}

