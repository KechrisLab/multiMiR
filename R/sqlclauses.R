

sql_org <- function(.table) { 
    where <- where_org(.table)
    as_mmsql(.where_list = as_where_list(where))
}

where_org <- function(.table) {
    no_target <- .table %in% tables_wo_target()
    as_where(.vars     = if (no_target) "m.org" else c("m.org", "t.org"),
             .connect  = "AND",
             .operator = "=")
}

sql_conserved <- function(.table, org, predicted.site) {
    where <- where_conserved(.table, org, predicted.site)
    as_mmsql(.where_list = as_where_list(where))
}

where_conserved <- function(.table, org, predicted.site) {

    miranda_cut    <- if (org == "mmu") 0.566 else 0.57
    targetscan_cut <- if (predicted.site == "conserved") "'N'" else "'Y'"
    pita_cut       <- 0.9

    vars      <- switch(.table, targetscan = "i.conserved_site", "i.conservation")
    operator  <- switch(predicted.site, conserved = ">=", "<")
    operator  <- switch(.table, targetscan = "=", operator)
    cut_value <- switch(.table, 
                        miranda    = miranda_cut,
                        targetscan = targetscan_cut,
                        pita       = pita_cut)
    operator  <- paste(operator, cut_value)

    has_conserved <- (.table %in% conserved_tables() & predicted.site != "all")

    as_where(.vars     = if (!has_conserved) NULL else vars,
             .operator = if (!has_conserved) NULL else operator)

}

sql_cutoff <- function(.table, org, predicted.site, predicted.cutoff.type,
                       predicted.cutoff) {

    score_var    <- get_score_var(.table)
    cutoff_name  <- create_cutoff_name(.table, org, predicted.site)
    score_cutoff <- cutoff_to_score(.table, cutoff_name, predicted.cutoff.type,
                                    predicted.cutoff)
    clause_cutoff(.table, score_var, score_cutoff)

}

clause_cutoff <- function(.table, score_var, score_cutoff) {

    operator <- switch(.table, miranda = "<=", pita = "<=", targetscan = "<=", ">")
    where    <- as_where(.vars = score_var, .operator = operator) #, .value = score_cutoff)
    # NOTE: Check what the .value should be here. wasn't gettting same value
    # from example3 in vignette
    
    as_mmsql(.where = as_where_list(where),
             .orderby = as_orderby(.vars = score_var, .order = "DESC"))

}

create_cutoff_name <- function(.table, org, predicted.site) {
    suffix <- switch(predicted.site, conserved = "c1", nonconserved = "c0", NULL)
    paste(c(.table, org, suffix), collapse = ".")
}

cutoff_to_score <- function(.table, cutoff_name, predicted.cutoff.type, predicted.cutoff) {

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
                            "number of records in table ", .table, ". All records ",
                            "will be queried.\n")
        if (predicted.cutoff < count_min) message(too_small)
        if (predicted.cutoff > tbl_count) message(too_large)

        adj_pred_cutoff <- max(min(tbl_count, predicted.cutoff), count_min)
        score_cutoff    <- cutoffs[[as.character(adj_pred_cutoff)]]
    }

    return(score_cutoff)

}

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


