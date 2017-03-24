
as_where_list <- function(...) {
    UseMethod("as_where_list")
}

as_where_list.default <- function(...) {
    wlist <- list(...)
    stopifnot(all(purrr::map_lgl(wlist[[1]], inherits, "where")))
    structure(..., class = "where_list")
}

as_where_list.where <- function(...) {
    wlist <- list(...)
    stopifnot(all(purrr::map_lgl(wlist, inherits, "where")))
    structure(wlist, class = "where_list")
}

# Constructor function 'where' class; i.e. sql where clauses
as_where <- function(.vars, .connect = NULL, .operator, .value = "%s") {
    structure(list(.vars     = .vars,
                   .connect  = .connect,
                   .operator = .operator,
                   .value    = .value),
              class = c("mmsql", "where"))
}

as_orderby <- function(.vars, .order) {
    stopifnot(.order %in% c("ASC", "DESC"))
    structure(list(.vars = .vars, .order = .order),
              class = c("mmsql", "orderby"))
}

as_mmsql <- function(.select = NULL, .from = NULL, .on = NULL, 
                     .where_list = NULL, .orderby = NULL) {

    if (!is.null(.where_list)) stopifnot(inherits(.where_list, "where_list"))
    if (!is.null(.orderby)) stopifnot(inherits(.orderby, "orderby"))

    structure(list(.select     = .select,
                   .from       = .from,
                   .on         = .on,
                   .where_list = .where_list,
                   .orderby    = .orderby),
              class = c("mmsql"))

}


################################################################################

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
    purrr::map(where_list, expand_where)
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

expand_sql <- function(x) {
# combine expand_x()'s
}

expand_select <- function(x) {
    paste("SELECT", paste(unlist(x), collapse = ", "))
}

expand_from <- function(x) {
    x <- merge_order(x)
    paste("FROM", paste(unlist(x), collapse = " INNER JOIN "))
}
expand_on <- function(x) {
    x <- merge_order(x)
    paste0("ON (", paste(unlist(x), collapse = " AND "), ")")
}
expand_wheres <- function(x) {
#NOTE: Left off here... its getting complicated, try to organize after figuring
    #out expanding wheres
    x <- 

}


# Reorder list to match the table merging order 
merge_order <- function(.list) {
    c(.list["mirna"],
      .list[c("validated", "predicted", "associations")],
      .list["target"])
}

build_mmquery <- function(.table, org, predicted.site, predicted.cutoff.type,
                          predicted.cutoff) { 

    features <- list(mirna        = sql_mirna(),
                     target       = sql_target(.table),
                     validated    = sql_validated(.table),
                     predicted    = sql_predicted(.table),
                     associations = sql_associations(.table),
                     org          = sql_org(.table),
                     conserved    = sql_conserved(.table, org, predicted.site),
                     cutoff       = sql_cutoff(.table, org, predicted.site,
                                               predicted.cutoff.type,
                                               predicted.cutoff))
    sql_parts <- transpose(features)
    sql_parts$.where_list <- 
        as_where_list(sql_parts$.where_list)
    sql_parts$.where_list %>% flatten %>% as_where_list %>% str
    paste(
          expand_select(sql_parts$.select),
          expand_from(sql_parts$.from),
          expand_on(sql_parts$.on)
          expand_where(sql_parts$.where_list)
          )

}



# just for saving code, not final form.
# expand_cutoff... <- function() {
#     paste(score_vars, operator, score_cutoff, "ORDER BY", score_vars, "DESC")
# }

