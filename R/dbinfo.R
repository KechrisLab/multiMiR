
#' Collect Information About the Web Server And Database of the multiMiR
#' Package
#' 
#' Functions for collecting and displaying information about the web server and
#' database of the multiMiR package.
#' 
#' \code{url} is a character string containing the URL of the multiMiR web
#' server.
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
#' \code{multimir_dbTables} returns the list of tables in the multiMiR
#' database.
#' 
#' @aliases multimir_dbCount multimir_dbInfo multimir_dbInfoVersions
#' multimir_dbSchema multimir_dbTables
#' @param url A character string naming the URL of the web server hosting the
#' multiMiR database.
#' @param schema.file A character string naming the file containing the
#' multiMiR database schema.
#' @return \code{multimir_dbCount}: a data frame with the count of records in
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
#' @author Yuanbin Ru \email{ruyuanbin@@gmail.com}
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
#' @export multimir_dbInfo
multimir_dbInfo <- function(url = getOption("multimir.url")) {
    qry <- "SELECT * FROM map_metadata"
    search.multimir(url = url, query = qry)
}

# To display database information on DB versions
multimir_dbInfoVersions <- function(url = getOption("multimir.url")) {
    qry <- paste("SELECT * FROM multimir_versions.version",
                 "WHERE public=1 ORDER BY version DESC")
    search.multimir(url = url, query = qry)
}

# To display database schema
multimir_dbSchema <- function(schema.file = getOption("multimir.schema.url")) {
    schema <- readLines(schema.file)
    cat(schema, sep = "\n")
}

# To show tables in the multimir database
multimir_dbTables <- function(url = getOption("multimir.db.tables")) {
    readLines(url)
}

