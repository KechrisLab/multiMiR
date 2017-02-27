
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
#' @param url A character string naming the URL of the web server hosting the
#' multiMiR database.
#' @param dbVer A character string containing the full version number for the
#' database version to use for for all package functions. The default will be
#' the most recent version.
#' @return
#' 
#' \code{multimir_dbInfoVersions}: a data frame with information about the
#' multiMiR database versions.
#' 
#' \code{multimir_switchDBVersion}: none (invisible \code{NULL}).
#' @author Yuanbin Ru \email{ruyuanbin@@gmail.com}
#' @keywords utilities database
#' @examples
#' 
#'   multimir_dbInfoVersions()
#'   multimir_switchDBVersion(dbVer="2.0.0")
#' 
#' @export multimir_switchDBVersion
multimir_switchDBVersion <- function(url = getOption("multimir.url"), 
                                     dbVer = getOption("multimir.db.version")) {
    # To switch DB version to search to the specified version if one matches
  tryCatch({
    query <- paste0("Select * from multimir_versions.version where version='",
                    dbVer, "' and public=1 order by version DESC")
    result <- RCurl::postForm(url, query = query, .cgifields = c("query"))
    result <- XML::readHTMLTable(result)

    if(as.numeric(as.character(result[[1]][[1]]))==0) {
      message("Version not set.\nVersion probably doesn't match an ",
              "available version.  Please use the full version as displayed in ",
              "multimir_dbInfoVersions().")
    } else {
      current <- result[[2]][1,]
      options(multimir.db.version  = as.character(current[[1]]))
      options(multimir.db.updated  = as.character(current[[2]]))
      options(multimir.db.name     = as.character(current[[4]]))
      options(multimir.db.tables   = paste0("http://multimir.ucdenver.edu/", 
                                            as.character(current[[7]])))
      options(multimir.url = url)
      options(multimir.schema.url  = paste0("http://multimir.ucdenver.edu/", 
                                            as.character(current[[5]])))
      options(multimir.cutoffs.url = paste0("http://multimir.ucdenver.edu/", 
                                            as.character(current[[3]])))
      options(multimir.error.msg   = "")
      cat(paste0("Now using database version: ", 
                  getOption("multimir.db.version")))
  }}, warning = function(war) {
      message(war)
  }, error = function(e) {
      message(e)
  }, finally = {})
  
}
