#' S3 Class constructors for objects defining SQL query components
#' and a collection of these parts (\code{mmsql_components}). 
#'
#' The collection object has a defined set of components that match the multiMiR
#' database and defined options in \code{get_multimir()}. Conceptually this is
#' split into two parts, the relatively straightforward SELECT, FROM, and ON
#' portion of the query and the more complex filtering and sorting operations:
#' WHERE and ORDER BY.  The latter have their own classes, the former are
#' resolved as strings (or character vectors) in the functions defining handling
#' of each sql table (\code{sql_} prefix).
#' 
#' @aliases as_mmsql_components, as_where_list, as_where, as_orderby,
#' is_where_list
#' @return 
#' \code{as_mmsql_components}: A collection of components that make up a
#' SQL query. 
#' \code{as_where_list}, \code{as_where}, \code{as_orderby}: Individual
#' components of a SQL query. 
#' @keywords internal
as_mmsql_components <- function(.select = NULL, .from = NULL, .on = NULL, 
                                .where_list = NULL, .orderby = NULL,
                                typeattr = NULL) {

    if (!is.null(.where_list)) stopifnot(inherits(.where_list, "where_list"))
    if (!is.null(.orderby)) stopifnot(inherits(.orderby, "orderby"))

    structure(list(.select     = .select,
                   .from       = .from,
                   .on         = .on,
                   .where_list = .where_list,
                   .orderby    = .orderby),
              class = c("mmsql_components"))

}

#' @rdname as_mmsql_components
#' @keywords internal
#' @importFrom purrr compact
#' @importFrom purrr map_lgl
as_where_list <- function(...) {
    wlist <- purrr::compact(list(...))
    stopifnot(all(purrr::map_lgl(wlist, inherits, "where")))
    structure(wlist, class = "where_list")
}

#' @rdname as_mmsql_components
#' @keywords internal
as_where <- function(.vars, .connect = NULL, .operator, .value = "%s") {
    structure(list(.vars     = .vars,
                   .connect  = .connect,
                   .operator = .operator,
                   .value    = .value),
              class = c("where"))
}

#' @rdname as_mmsql_components
#' @keywords internal
is_where <- function(x) inherits(x, "where")

#' @rdname as_mmsql_components
#' @keywords internal
as_orderby <- function(.vars, .order) {
    stopifnot(.order %in% c("ASC", "DESC"))
    structure(list(.vars = .vars, .order = .order),
              class = c("orderby"))
}


