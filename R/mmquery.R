
#' Constructor and methods for object returned by \code{get.multimir()}.
#' 
#' This package's primary user-facing object. Contains the SQL statement and the
#' returned data query, as well as summary and link objects depending on
#' specified options 
#' 
#' @keywords internal
as_mmquery <- function(.data, queries, summary = NULL, links = NULL) {
    structure(list(validated,
                   predicted,
                   disease.drug,
                   query = query, #list(query, sql_parts),
                   summary,
                   links),
              class = c("mmquery"),
              types = c(),
              tables = c())
}
# Or
as_mmquery <- function(){

    sqlvalidated   <- query[[names(query) %in% validated_tables()]]
    sqlpredicted   <- query[[names(query) %in% predicted_tables()]]
    sqldiseasedrug <- query[[names(query) %in% diseasedrug_tables()]]
    attr(validated, "query")    <- sqlvalidated
    attr(predicted, "query")    <- sqlpredicted
    attr(diseasedrug, "query")  <- sqldiseasedrug

    structure(list(validated = validated,
                   predicted = predicted,
                   disease.drug = diseasedrug,
                   summary,
                   links),
              class = c("mmquery"))
}

print.mmquery <- function() {

}

summary.mmquery <- function() {

}



