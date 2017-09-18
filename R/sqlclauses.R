#' Functions defining the WHERE clauses.
#'
#' Functions defining filtering by organism (org), disease/drug, conserved, and
#' cutoff. Filtering by mirna and target are defined within their \code{sql_...}
#' functions.
#' 
#' @aliases sql_org, where_org, where_diseasedrug, where_conserved,
#' where_cutoff create_cutoff_name cutoff_to_score
#' @return The \code{WHERE} portion of a SQL query
#' @keywords tables types predicted validated diseasedrug disease drug
#' 
#' @keywords internal
sql_org <- function(.table, org) { 
    where <- where_org(.table, org)
    as_mmsql_components(.where_list = as_where_list(where))
}

#' @rdname sql_org
#' @keywords internal
where_org <- function(.table, org) {
    no_target <- .table %in% tables_wo_target()
    as_where(.vars     = if (no_target) "m.org" else c("m.org", "t.org"),
             .connect  = "AND",
             .operator = "=",
             .value    = quote_wrap(org))
}


#' @rdname sql_org
#' @keywords internal
where_diseasedrug <- function(.table, disease.drug) {
    wherevar <- switch(.table,
                       pharmaco_mir = "i.drug",
                       mir2disease = "i.disease", 
                       phenomir = c("i.disease", "i.disease_class"))
    if (is.null(disease.drug)) {
        NULL
    } else {
        as_where(.vars = wherevar, 
                 .connect  = "OR", 
                 .operator = "IN",
                 .value = disease.drug)
    }

}


#' @rdname sql_org
#' @keywords internal
where_conserved <- function(.table, org, predicted.site) {

    miranda_cut    <- switch(org, mmu = 0.566, 0.57)
    targetscan_cut <- if (predicted.site == "conserved") "'Y'" else "'N'"
    pita_cut       <- 0.9

    vars      <- switch(.table, targetscan = "i.conserved_site",
                        "i.conservation")
    operator  <- switch(predicted.site, conserved = ">=", "<")
    operator  <- switch(.table, targetscan = "=", operator)
    cut_value <- switch(.table, 
                        miranda    = miranda_cut,
                        targetscan = targetscan_cut,
                        pita       = pita_cut)

    has_conserved <- (.table %in% conserved_tables() & predicted.site != "all")

    if (!has_conserved) 
        NULL
    else {
        as_where(.vars = vars, .operator = operator, .value = cut_value)
    }

}


#' @rdname sql_org
#' @keywords internal
where_cutoff <- function(.table, score_var, score_cutoff) {

    # NOTE: Check what the .value should be here. wasn't gettting same value
    # from example3 in vignette
    operator <- switch(.table, miranda = "<=", pita = "<=", targetscan = "<=",
                       ">")
    as_where(.vars = score_var, .operator = operator, .value = score_cutoff)

}


#' @rdname sql_org
#' @keywords internal
create_cutoff_name <- function(.table, org, predicted.site) {
    if (.table %in% conserved_tables()) {
        suffix <- switch(predicted.site, conserved = "c1", nonconserved = "c0",
                         NULL)
    } else {
        suffix <- NULL
    }
    paste(c(.table, org, suffix), collapse = ".")
}


#' @rdname sql_org
#' @keywords internal
cutoff_to_score <- function(.table, cutoff_name, predicted.cutoff.type,
                            predicted.cutoff) {

    scipen.orig <- getOption("scipen")
    options(scipen = 99)
    on.exit(options(scipen = scipen.orig))
    cutoffs   <- get.multimir.cutoffs()[[cutoff_name]]

    # get dataset-specific score cutoff
    if (predicted.cutoff.type == "p") {

        score_cutoff <- cutoffs[[paste0(predicted.cutoff, "%")]]

    } else if (predicted.cutoff.type == "n") {

        tbl_count <- cutoffs[["count"]]

        count_min        <- 10000
        too_small <- paste("Number predicted cutoff (predicted.cutoff)",
                           predicted.cutoff, "may be too small. A cutoff of",
                           "10000 will be used instead.\n")
        too_large <- paste0("Number predicted cutoff (predicted.cutoff) ",
                            predicted.cutoff, " is larger than the total ",
                            "number of records in table ", .table, 
                            ". All records will be queried.\n")
        if (predicted.cutoff < count_min) message(too_small)
        if (predicted.cutoff > tbl_count) message(too_large)

        adj_pred_cutoff <- max(min(tbl_count, predicted.cutoff), count_min)
        adj_pred_cutoff <- as.integer(as.integer(adj_pred_cutoff / 10000) *
                                      10000)
        score_cutoff    <- cutoffs[[as.character(adj_pred_cutoff)]]
    }

    return(score_cutoff)

}



