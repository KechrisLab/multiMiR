
this_table  <- 
    #"diana_microt" %>% 
    #"mir2disease" %>% 
    "targetscan" %>% sql_features
select_list <- this_table %>% transpose %>% .$.select
from_list   <- this_table %>% transpose %>% .$.from
on_list     <- this_table %>% transpose %>% .$.on
where_list  <- this_table %>% transpose %>% .$.where

paste(expand_select(select_list),
      expand_from(from_list),
      expand_on(on_list))

combine_wheres(where_list)

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

sql_features <- function(.table) { 

    target <- sql_target(.table)

    list(mirna        = sql_mirna(),
         target       = target,
         validated    = sql_validated(.table),
         predicted    = sql_predicted(.table),
         associations = sql_associations(.table),
         org          = sql_org(target$.where))

    }

sql_validated <- function(.table) {

    these_tables <- c("mirecords", "mirtarbase", "tarbase")
    this_type <- .table %in% these_tables

    .select <- if (this_type) c("i.experiment, i.support_type, i.pubmed_id")  else NULL
    .from   <- if (this_type) "%s AS i" else NULL
    .on     <- NULL
    .where  <- list(vars = NULL, connect  = NULL, operator = NULL, value = NULL)

    return(list(.select = .select, .from = .from, .on = .on, .where = .where))
}

sql_predicted <- function(.table) {

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
    .where  <- list(vars = NULL, connect  = NULL, operator = NULL, value = NULL)

    return(list(.select = .select, .from = .from, .on = .on, .where = .where))
}

sql_associations <- function(.table) {
    
    these_tables <- c("mir2disease", "pharmaco_mir", "phenomir")
    this_type <- .table %in% these_tables

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

    .select <- if (!this_type) NULL else .select
    .from   <- if (!this_type) NULL else sprintf("%s AS i", .table)
    .on     <- NULL
    .where  <- TRUE
    .where  <- list(vars = wherevar, connect  = "OR", operator = "IN", value = "%s")

    return(list(.select = .select, .from = .from, .on = .on, .where = .where))
}

sql_mirna <- function() {
    .select <- c("m.mature_mirna_acc, m.mature_mirna_id")
    .from   <- "mirna AS m"
    .on     <- "m.mature_mirna_uid = i.mature_mirna_uid"
    .where  <- list(vars     = c("m.mature_mirna_acc", "m.mature_mirna_id"),
                    connect  = "OR",
                    operator = "IN",
                    value    = "i.mature_mirna_uid")
    return(list(.select = .select, .from = .from, .on = .on, .where = .where))
}

sql_org <- function(target_where) {

    no_target <- if (is.null(target_where$vars)) TRUE else FALSE
    .select <- NULL
    .from   <- NULL
    .on     <- NULL
    .where  <- list(vars     = (if (no_target) "m.org" else c("m.org", "t.org")),
                    connect  = "AND",
                    operator = "=",
                    value    = "%s")

    return(list(.select = .select, .from = .from, .on = .on, .where = .where))
}

sql_target <- function(.table) {

    tables_wo_target <- c("mir2disease", "phenomir")
    no_target        <- .table %in% tables_wo_target

    vars    <- c("t.target_symbol", "t.target_entrez", "t.target_ensembl")
    .select <- paste(paste0("%1$s", vars), collapse = ", ")
    na_txt  <- "'NA' AS "

    .select <- if (no_target) sprintf(.select, na_txt) else sprintf(.select, "")
    .from   <- if (no_target) NULL else "target AS t"
    .on     <- if (no_target) NULL else "i.target_uid = t.target_uid"
    .where  <- list(vars     = if (no_target) NULL else vars,
                    connect  = "OR",
                    operator = "IN",
                    value    = "%s")

    return(list(.select = .select, .from = .from, .on = .on, .where = .where))

}

expand_select <- function(select) {
    paste("SELECT", paste(unlist(select), collapse = ", "))
}
expand_from <- function(from) {
    from <- merge_order(from)
    paste("FROM", paste(unlist(from), collapse = ", "))
}
expand_on <- function(on) {
    on <- merge_order(on)
    paste0("ON (", paste(unlist(on), collapse = " AND "), ")")
}



where_conserved <- function(.table, predicted.site, org) {

    conserved_tables <- c("miranda", "pita", "targetscan")
    has_conserved    <- (.table %in% conserved_tables & predicted.site != "all")

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

    .where <- list(vars     = if (!has_conserved) NULL else vars,
                   connect  = NULL,
                   operator = if (!has_conserved) NULL else operator,
                   value    = "%s")

    return(list(.select = NULL, .from = NULL, .on = NULL, .where = .where))

}

where_cutoff <- function(.table, score_vars, score_cutoff) {

    operator <- switch(.table,
                       miranda = "<=",
                       pita = "<=",
                       targetscan = "<=",
                       ">")
    qry      <- paste(score_vars, operator, score_cutoff, "ORDER BY", score_vars, "DESC")

}

cutoff_name <- function(.table, predicted.site, org) {
    suffix <- switch(predicted.site, conserved = "c1", nonconserved = "c0", NULL)
    paste(c(table, org, suffix), collapse = ".")
}

cutoff_to_score <- function(cutoff_name, predicted.cutoff.type, predicted.cutoff) {

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
                            "number of records in table ", table, ". All records ",
                            "will be queried.\n")
        if (predicted.cutoff < count_min) message(too_small)
        if (predicted.cutoff > tbl_count) message(too_large)

        adj_pred_cutoff <- max(min(tbl_count, predicted.cutoff), count_min)
        score_cutoff    <- cutoffs[[as.character(adj_pred_cutoff)]]
    }

    return(score_cutoff)

}


# Expanding the simple wheres is a multi-step process
# This expands a single where list (i.e. org()$.where)
expand_where <- function(where) {
    # note, each sql feature has its own where clause
    sep <- pad(where$connect)
    if (is.null(where$vars)) {
        return(NULL)
    } else {
        where <- purrr::map_chr(where$vars, paste, where$operator, where$value)
        paste0("(", paste(where, collapse = sep), ")")
    }
}

# Maps expand_where to all single wheres in where_list 
# note: this could combine with where_args to expand full where clause (of
# simple ones at least -- excluding conserved and cutoff
combine_wheres <- function(string_list) {
    #.lst <- 
    map(where_list, expand_where)
    #unlist(.lst)
}

# For select and adding user input args to the simple where statements
where_args <- function(mirna, target, disease.drug, org) {
    list(mirna        = mirna,
         target       = target,
         validated    = NULL,
         predicted    = NULL,
         associations = disease.drug,
         org          = org)
}


# Reorder list to match the table merging order 
merge_order <- function(.list) {
    c(.list["mirna"],
      .list[c("validated", "predicted", "associations")],
      .list["target"])
}

# pad single space on each side
pad <- function(x) paste0(" ", x, " ") 
