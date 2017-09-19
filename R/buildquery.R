#' Constructors for parts of SQL queries
#' Expand_query converts a \code{mmyquery} object to a SQL query string
#' 
#' @return A complete SQL statement and related information.
#' @keywords internal
#' @importFrom purrr map
#' @importFrom purrr transpose
#' @importFrom purrr compact
#' @importFrom purrr map_chr
#' @importFrom purrr flatten
build_mmsql <- function(.table, org, 
                        mirna                 = NULL,
                        target                = NULL,
                        disease.drug          = NULL, 
                        predicted.site        = NULL,
                        predicted.cutoff.type = NULL,
                        predicted.cutoff      = NULL,
                        limit                 = NULL) {

    components <- list(mirna        = sql_mirna(mirna),
                       target       = sql_target(.table, target),
                       validated    = sql_validated(.table),
                       predicted    = sql_predicted(.table, org, predicted.site,
                                                    predicted.cutoff.type,
                                                    predicted.cutoff),
                       diseasedrug  = sql_diseasedrug(.table, disease.drug),
                       org          = sql_org(.table, org))

    sql_parts       <- purrr::transpose(components)
    sql_parts$.limit <- limit
    sql_parts_trim  <- purrr::map(sql_parts, purrr::compact)
    table_type      <- reverse_table_lookup(.table)

    
    list(query = expand_query(sql_parts_trim),
         #sqlparts = sql_parts_trim,
         table = .table,
         type  = table_type)

}

#' @rdname build_mmsql
#' @keywords internal
expand_query <- function(x) {
    paste(expand_select(x$.select),
          expand_from(x$.from),
          expand_on(x$.on),
          expand_where_list(x$.where_list),
          expand_orderby(x$.orderby),
          expand_limit(x$.limit))
}

#' @rdname build_mmsql
#' @keywords internal
expand_select <- function(x) {
    paste("SELECT", paste(unlist(x), collapse = ", "))
}

#' @rdname build_mmsql
#' @keywords internal
expand_from <- function(x) {
    x <- merge_order(x)
    paste("FROM", paste(unlist(x), collapse = " INNER JOIN "))
}

#' @rdname build_mmsql
#' @keywords internal
expand_on <- function(x) {
    x <- merge_order(x)
    paste0("ON (", paste(unlist(x), collapse = " AND "), ")")
}

#' @rdname build_mmsql
#' @keywords internal
expand_where_list <- function(x) {
    paste("WHERE", paste(unlist(purrr::map(purrr::flatten(x), expand_where)),
                         collapse = " AND "))

}

#' @rdname build_mmsql
#' @keywords internal
expand_where <- function(x) {
    sep <- pad(x$.connect)
    value <- switch(x$.operator,
                    IN = parens_quote(x$.value),
                    x$.value)
    x <- purrr::map_chr(x$.vars, paste, x$.operator, value)
    parens_wrap(paste(x, collapse = sep))
}

#' @rdname build_mmsql
#' @keywords internal
expand_orderby <- function(x) {
    if (!is.null(x) & length(x) > 0) {
        x <- purrr::compact(x)
        stopifnot(length(x) == 1, inherits(x[[1]], "orderby"))
        x <- purrr::flatten(x)
        paste("ORDER BY", x$.vars, x$.order)
    }
}

#' @rdname build_mmsql
#' @keywords internal
expand_limit <- function(x) {
    if (!is.null(x)) paste("LIMIT", x)
}

#' @rdname build_mmsql
#' @keywords internal
merge_order <- function(.list) {
    # Reorder list to match the table merging order 
    c(.list["mirna"],
      .list[c("validated", "predicted", "diseasedrug")],
      .list["target"])
}

