
#' Collect Information About the Web Server And Database of the multiMiR
#' Package
#' 
#' Functions for collecting and displaying information about the web server and
#' database of the multiMiR package.
#' 
#' \code{multimir.url} is a global option containing the URL of the multiMiR web
#' server. Set using \code{options("multimir.url" = ...)}
#' 
#' \code{multimir_dbCount} returns counts of records in the tables in the
#' multiMiR database. Each table contains data from an external miRNA/target
#' database.
#' 
#' \code{multimir_dbInfo} returns other information about the multiMiR
#' database. This includes information of external miRNA/target databases in
#' multiMiR.
#' 
#' \code{multimir_dbInfoVersions} returns other information about the multiMiR
#' database versions available. This provides a list of available options if
#' switching to previous version is desired.
#' 
#' \code{multimir_dbSchema} prints the schema definition of the multiMiR
#' database.
#' 
#' \code{multimir_dbTables} returns the vector of tables in the multiMiR
#' database and saves it to the global option \code{multimir.tables.list}. This
#' function is automatically run when get_multimir is called if the
#' \code{multimir.tables.list} is NULL.
#' 
#' @aliases multimir_dbCount multimir_dbInfo multimir_dbInfoVersions
#' multimir_dbSchema multimir_dbTables
#' @param url Deprecated. Use global option \code{multimir.url} instead. 
#' @param schema.file Deprecated. Option exists as \code{multimir.schema},
#' but it should not need to be set directly.
#' @return 
#' \code{multimir_dbCount}: a data frame with the count of records in
#' each of the tables in the multiMiR database.
#' 
#' \code{multimir_dbInfo}: a data frame with information about the multiMiR
#' database.
#' 
#' \code{multimir_dbInfoVersions}: a data frame with information about the
#' multiMiR database versions.
#' 
#' \code{multimir_dbSchema}: none (invisible \code{NULL}).
#' 
#' \code{multimir_dbTables}: a data frame with table names in the multiMiR
#' database.
#' @keywords utilities database
#' @examples
#' 
#'   this_url <- getOption("multimir.url")
#'   this_url
#'   options(multimir.url = this_url)
#'   
#'   db_ver <- multimir_dbInfoVersions()
#'   
#'   db_count <- multimir_dbCount()
#' 
#'   db_info <- multimir_dbInfo()
#' 
#'   multimir_dbSchema()
#' 
#'   db_tables <- multimir_dbTables()
#' 
#' @export 
multimir_dbInfo <- function(url = NULL) {

    if (!is.null(url)) deprecate_arg("url")
    qry <- "SELECT * FROM map_metadata"
    search_multimir(query = qry)

}


#' @rdname multimir_dbInfo
#' @export
multimir_dbInfoVersions <- function(url = NULL) {

    if (!is.null(url)) deprecate_arg("url")
    qry <- paste("SELECT * FROM multimir_versions.version",
                 "WHERE public=1 ORDER BY version DESC")
    submit_request(query = qry, .cgifields = "query")
}

#' @rdname multimir_dbInfo
#' @export
multimir_dbSchema <- function(schema.file = NULL) {

    if (!is.null(schema.file)) deprecate_arg("schema.file")
    schema <- readLines(full_url("multimir.schema"))
    cat(schema, sep = "\n")

}

#' @rdname multimir_dbInfo
#' @export
multimir_dbTables <- function(url = NULL) {

    if (!is.null(url)) deprecate_arg("db.tables")
    tbls <- as.character(readLines(full_url("multimir.db.tables")))
    options("multimir.tables.list" = tbls)
    return(tbls)

}

#' @rdname multimir_dbInfo
#' @export
multimir_dbCount <- function(url = NULL) {

    if (!is.null(url)) deprecate_arg("url")

    qry <- "SELECT * FROM map_counts"
    res <- search_multimir(query = qry)

    for (i in 2:ncol(res)) {
        res[, i] <- as.numeric(as.character(res[, i]))
    }

    return(res)

}

