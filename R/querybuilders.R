#' Prepare query for validated target table
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get.multimir}.
#'
#' @param table One of the three validated tables
#' @keywords internal
query_validated <- function(table = c("mirecords", "mirtarbase", "tarbase"), 
                            mirna = NULL, target = NULL, org, mirna.table,
                            target.table, ...) {

	table <- match.arg(table)

    paste("SELECT m.mature_mirna_acc, m.mature_mirna_id,",
          "t.target_symbol, t.target_entrez, t.target_ensembl,",
          "i.experiment, i.support_type, i.pubmed_id",
          "FROM", mirna.table, "AS m INNER JOIN", table, 
          "AS i INNER JOIN", target.table, 
          "AS t ON (m.mature_mirna_uid=i.mature_mirna_uid AND",
          "i.target_uid=t.target_uid) WHERE",
          subquery_mirnatarget(mirna, target), 
          "AND",
          subquery_org(table, org))

}


#' Prepare query for predicted target table
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get.multimir}.
#'
query_predicted <- function(table = c("diana_microt", "elmmo", "microcosm",
                                      "miranda", "mirdb", "pictar", "pita",
                                      "targetscan"), 
                            mirna = NULL, target = NULL, org, mirna.table,
                            target.table, predicted.cutoff,
                            predicted.cutoff.type, 
                            predicted.site = c("conserved", "nonconserved", "all"),
                            ...) {
	table <- match.arg(table)
	predicted.site <- match.arg(predicted.site)

    # NOTE: Should this be a message saying org 'rno' doesn't exist for these tables?
    if (org == "rno" & table %in% c("diana_microt", "pictar", "pita",
                                    "targetscan")) return(NULL)

    # Set name of score variable by table
    score_vars <- switch(table,
                         diana_microt = "i.miTG_score",
                         elmmo        = "i.p ",
                         microcosm    = "i.score", # NOTE: didn't have the 'AS score' text
                         mirdb        = "i.score", #       didn't have the 'AS score' text
                         pictar       = "i.score", #       didn't have the 'AS score' text
                         miranda      = "i.mirsvr_score",
                         pita         = "i.ddG",
                         targetscan   = "i.context_plus_score")

    # Create subquery_conserved and 'name' variable 
    conserved_list <- subquery_conserved(predicted.site = predicted.site, 
                                        table = table, org = org)


    # Basic SQL structure:
    # predicted: X(choose score var) + Xbase + X(mirna and/or target) + Xorg + #
    #            Xconserved(for 3 of 8) + Xcutoff
	query_list <- 
        list(base = paste("SELECT m.mature_mirna_acc, m.mature_mirna_id,",
                          "t.target_symbol, t.target_entrez, t.target_ensembl,",
                          score_vars,
                          "AS score FROM", mirna.table, "AS m INNER JOIN", table, 
                          "AS i INNER JOIN", target.table, 
                          "AS t ON (m.mature_mirna_uid=i.mature_mirna_uid",
                          "AND i.target_uid=t.target_uid) WHERE"),
             mirna_target = subquery_mirnatarget(mirna, target),
             org          = subquery_org(table = table, org = org),
             conserved    = conserved_list$query,
             cutoff       = subquery_cutoff(table, predicted.cutoff.type, 
                                            predicted.cutoff,
                                            conserve_list$name, score_vars))

    query_list <- query_list[do.call(c, lapply(query_list, function(x) !is.na(x)))]
    qry        <- paste(query_list[[1]], paste(do.call(c, query_list[-1]),
                                               collapse = " AND ")) 

	return(qry)

}


#' Prepare query for miRNA-disease/drug table
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get.multimir}.
#'
query_disease <- function(table = c("mir2disease", "pharmaco_mir", "phenomir"), 
						  mirna = NULL, target = NULL, org, mirna.table, target.table,
						  disease.drug, ...) {

	table <- match.arg(table)
    if (table %in% c("mir2disease", "phenomir") & 
        is.null(mirna) & is.null(disease.drug)) return(NULL)

    query_base <- 
        switch(table,
               mir2disease = paste("SELECT m.mature_mirna_acc, m.mature_mirna_id, 'NA' AS",
                                   "target_symbol, 'NA' AS target_entrez, 'NA' AS",
                                   "target_ensembl, i.disease AS disease_drug,",
                                   "CONCAT_WS('. ', i.year, i.title)",
                                   "AS paper_pubmedID FROM", 
                                   mirna.table, "AS m INNER JOIN", table, 
                                   "AS i ON (m.mature_mirna_uid=i.mature_mirna_uid) WHERE"),
               pharmaco_mir = paste("SELECT m.mature_mirna_acc, m.mature_mirna_id,",
                                    "t.target_symbol, t.target_entrez, t.target_ensembl,",
                                    "i.drug AS disease_drug, i.pubmed_id AS paper_pubmedID",
                                    "FROM", mirna.table, "AS m INNER JOIN", table, 
                                    "AS i INNER JOIN", target.table, 
                                    "AS t ON (m.mature_mirna_uid=i.mature_mirna_uid AND",
                                    "i.target_uid=t.target_uid) WHERE"),
               phenomir     = paste("SELECT m.mature_mirna_acc, m.mature_mirna_id, 'NA' AS",
                                    "target_symbol, 'NA' AS target_entrez, 'NA' AS",
                                    "target_ensembl, i.disease AS disease_drug, i.pubmed_id AS",
                                    "paper_pubmedID FROM", mirna.table, "AS m INNER JOIN", table, 
                                    "AS i ON (m.mature_mirna_uid=i.mature_mirna_uid) WHERE"))
    query_diseasedrug <- switch(table, 
                              mir2disease  = paste("i.disease IN", disease.drug),
                              pharmaco_mir = paste("i.drug IN", disease.drug),
                              phenomir     = paste("(i.disease IN", disease.drug, "OR",
                                                   "i.disease_class IN", disease.drug, ")"))
    query_diseasedrug <- ifelse(is.null(disease.drug), "", query_diseasedrug)

    query_org    <- subquery_org(table = table, org = org)
    query_target <- ifelse(table != "pharmaco_mir", "",
                           subquery_mirnatarget(target = target))
    query_mirna  <- subquery_mirnatarget(mirna = mirna)
    
    # Basic SQL structure is:
    #   mir2disease: base + disease.drug + org + mirna
    #   pharmaco_mir: base + target + disease.drug + org + mirna
    #   phenomir: base +  disease.drug + org + mirna
    query_list <- list(query_base, query_target, query_diseasedrug, query_org,
                       query_mirna)
    query_list <- query_list[do.call(c, lapply(query_list, function(x) x != ""))]
    rtn        <- paste(query_list[[1]], paste(do.call(c, query_list[-1]),
                                               collapse = " AND ")) 

	return(rtn)
}



#' Internal function for building mirna/target portion of query
#'
#' @param mirna see ?get.multimir
#' @param target see ?get.multimir
#' @keywords internal
subquery_mirnatarget <- function(mirna = NULL, target = NULL) {

    if(is.null(mirna) & is.null(target)) {
        stop("Args mirna and target are NULL. At least one must be specified ",
             "for the selected table")
    }

    query_list <- 
        list(mirna  = paste("( m.mature_mirna_acc IN", mirna, 
                            "OR m.mature_mirna_id IN", mirna, ")"),
             target = paste("( t.target_symbol IN", target, 
                            "OR t.target_entrez IN", target, 
                            "OR t.target_ensembl IN", target, ")"))
    query_logic <- 
        list(mirna          = !is.null(mirna),
             target         = !is.null(target))

	query_list <- query_list[do.call(c, query_logic)]
	do.call(paste, c(query_list, sep = " AND "))

}

#' Internal function for building org portion of query
#' 
#' @param table One of the validate, predict, or disease tables
#' @param org One of the 3 organisms (hsa, mmu, rno)
#' @keywords internal
subquery_org <- function(table, org) {

    stopifnot(!is.null(org), !is.null(table))

    if (table %in% c("tarbase", "mir2disease", "phenomir")) {
        paste("m.org =", quote_wrap(org))
    } else {
		paste("m.org =", quote_wrap(org), "AND t.org =", quote_wrap(org))
    }

}

#' Internal function for building conserved (predicted.site) portion of query
#' 
#' Process conserved/nonconserved sites in miranda, pita & targetscan
#' (other tables don't have conserved/nonconserved annotation)
#'   - 2 objects created: 'name' and 'query'
#'   - subquery_conserved structure: conserved_var + operator + cut_value
#'
#' @param predicted.site see ?get.multimir
#' @param table see ?get.multimir
#' @param org see ?get.multimir
#' @keywords internal
subquery_conserved <- function(predicted.site, table, org) {

    if (table %in% c("miranda", "pita", "targetscan")) {
        conserv_args <- 
            list(suffix    = switch(predicted.site,
                                    conserved    = "c1",
                                    nonconserved = "c0"),
                 variable  = switch(table,
                                    pita = "i.conservation",
                                    miranda = "i.conservation",
                                    targetscan = "i.conserved_site"),
                 operator  = switch(table,
                                    pita       = ifelse(predicted.site == "conserved", ">=", "<"),
                                    miranda    = ifelse(predicted.site == "conserved", ">=", "<"),
                                    targetscan = "="),
                 cut_value = switch(table,
                                    pita       = "0.9",
                                    miranda    = ifelse(org == "mmu", "0.566", "0.57"),
                                    targetscan = "'N'"))

        name <- paste(table, org, conserv_args$suffix, sep = ".")
        qry  <- paste(do.call(c, conserv_args[-1]), collapse = " ")
    } else {
        name <- paste(table, org, sep = ".")
        qry  <- ""
    }

    return(list(`name` = name, `query` = qry))

}

#' Internal function for creating the portion of the query for filtering based
#' on score
#' 
#' Depends on query_predicted() for input arguments.
#'
#' @param table
#' @param predicted.cutoff.type
#' @param predicted.cutoff
#' @param name
#' @param score_vars
#' @keywords internal
subquery_cutoff <- function(table, predicted.cutoff.type, predicted.cutoff,
                            name, score_vars) {

	cutoffs          <- get.multimir.cutoffs()
    tbl_count        <- cutoffs[[name]][["count"]]
    count_min        <- 10000
    count_max        <- 300000
    cutoff_too_small <- paste("Number predicted cutoff (predicted.cutoff)",
                              predicted.cutoff, "may be too small. A cutoff",
                              "of 10000 will be used instead.\n")
    cutoff_too_large <- paste0("Number predicted cutoff (predicted.cutoff) ",
                               predicted.cutoff, " is larger than the total ",
                               "number of records in table ", table, 
                               ". All records will be queried.\n")

	# get dataset-specific score cutoff
	if (predicted.cutoff.type == "p") {
        score.cutoff <- cutoffs[[name]][[paste0(predicted.cutoff, "%")]]
	} else if (predicted.cutoff.type == "n") {
        if (predicted.cutoff < count_min) message(cutoff_too_small)
        if (predicted.cutoff > tbl_count) message(cutoff_too_large)

        cut_pred.cutoff <- cut(predicted.cutoff, labels = FALSE,
                               breaks = c(0, count_min, tbl_count))
        adj.pred.cutoff <- switch(as.character(cut_pred.cutoff), 
                                  '1' = count_min, '2' = predicted.cutoff, 
                                  '3' = tbl_count, NA)
        adj.pred.cutoff <- ifelse(tbl_count < adj.pred.cutoff & 
                                  tbl_count >= count_min, 
                                  tbl_count, adj.pred.cutoff)
        score.cutoff    <- ifelse(is.na(adj.pred.cutoff), NA,
                                  cutoffs[[name]][[as.character(adj.pred.cutoff)]])
    }

    operator <- ifelse(table %in% c("miranda", "pita", "targetscan"), "<=", ">=")
    qry      <- ifelse(is.na(score.cutoff), NA, 
                       paste(score_vars, ">=", score.cutoff, "ORDER BY",
                             score_vars, "DESC"))
	return(qry)

}

