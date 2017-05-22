
#' Manage Database Version to use
#' 
#' Functions for managing the database version used to complete requests on the
#' web server.
#' 
#' \code{url} is a character string containing the URL of the multiMiR web
#' server. Optional as it is set when the package is loaded.
#' 
#' \code{multimir_dbInfoVersions} returns other information about the multiMiR
#' database versions available. This provides a list of available options if
#' switching to previous version is desired.
#' 
#' \code{multimir_switchDBVersion} returns other information about the multiMiR
#' database versions available. This provides a list of available options if
#' switching to previous version is desired.
#' 
#' @param url Deprecated. Use global option \code{multimir.url} instead. 
#' @param db_version A character string containing the full version number for
#' the database version to use for for all package functions. The default will
#' be the most recent version.
#' @return
#' 
#' \code{multimir_dbInfoVersions}: a data frame with information about the
#' multiMiR database versions.
#' 
#' \code{multimir_switchDBVersion}: none (invisible \code{NULL}).
#' @keywords utilities database
#' @examples
#' 
#'   multimir_dbInfoVersions()
#'   multimir_switchDBVersion(db_version="2.0.0")
#' 
#' @export multimir_switchDBVersion
multimir_switchDBVersion <- function(db_version, url = NULL) {

    if (!is.null(url)) deprecate_arg("url")

    old_vers   <- getOption("multimir.db.version")
    vers_table <- multimir_dbInfoVersions()
    new_vers   <- subset(vers_table, vers_table$VERSION == db_version)

    if (nrow(new_vers) == 0) stop("DB version ", db_version, " not found. ",
                                  "Please use the full version as displayed ",
                                  "in multimir_dbInfoVersions().")
    if (nrow(new_vers) > 1)  stop("Multiple matching DB versions.  Consider ",
                                  "submitting an issue on Github repository.")

    set_dbversion(dbversion_row = new_vers, overwrite = TRUE)
    message(paste0("Now using database version: ",
                   getOption("multimir.db.version")))

}



# Internal function
# @param overwrite Overwite existing options (TRUE/FALSE).
# @param dbversion_row A single row data frame from \code{queryDBVersions()}
set_dbversion <- function(dbversion_row, overwrite = FALSE) {

    # Get current global options
    op <- options()

    # Convert all cols to character (using [] is a trick via SO and Hadley that
    # allows the dataframe to keep its 'data.frame' class)
    dbversion_row[] <- lapply(dbversion_row, as.character)

    op.multimir.vers <-
        list(multimir.db.version = dbversion_row$VERSION,
             multimir.db.updated = dbversion_row$UPDATED,
             multimir.db.name    = dbversion_row$DBNAME,
             multimir.db.tables  = dbversion_row$TABLES,
             multimir.schema     = dbversion_row$SCHEMA,
             multimir.cutoffs    = dbversion_row$RDA)

    # Set option values
    toset <- ifelse(overwrite, 
                    rep(TRUE, length(op.multimir.vers)),
                    !(names(op.multimir.vers) %in% names(op)))
    if (any(toset)) options(op.multimir.vers[toset])

}


# Internal function for combining base url (multimir.url) with the url paths for
# various requests.
# @param pkg_option One of the package options corresponding to a url path
# needing appending to the base url 
# @return the full url to the object path passed to the function 
full_url <- function(pkg_option = c("multimir.queries",
                                    "multimir.db.tables",
                                    "multimir.schema",
                                    "multimir.cutoffs")) {

    pkg_option <- match.arg(pkg_option)
    paste0(getOption("multimir.url"), getOption(pkg_option))

}






