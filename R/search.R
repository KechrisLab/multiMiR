
#' Search the multiMiR Database Given a MySQL Query
#' 
#' This is an internal function called by some of the other functions in
#' multiMiR. Given a MySQL query, it searches and retrieves result from the
#' multiMiR database on the multiMiR web server. To use \code{search.multimir}
#' directly, users will need to be familiar with MySQL and multiMiR table
#' structures. Users are advised to use \code{get.multimir} instead.
#' 
#' 
#' @param url a character string for the URL of the multiMiR web server.  The
#' default is getOption("multimir.url")
#' ("http://multimir.ucdenver.edu/cgi-bin/multimir.pl").
#' @param query a character string for the MySQL query.
#' @return \code{search.multimir} returns a data frame containing results from
#' the multiMiR web server.
#' @author Yuanbin Ru \email{ruyuanbin@@gmail.com}
#' @keywords utilities database
#' @examples
#' 
#'   ## show all tables in the multiMiR database
#'   tables <- search.multimir(query="show tables")
#' 
#'   ## show the structure of table diana_microt
#'   microt <- search.multimir(query="describe diana_microt")
#' 
#'   ## search for validated target genes of hsa-miR-18a-3p in miRecords
#'   qry <- paste("SELECT m.mature_mirna_acc, m.mature_mirna_id,",
#'                "    t.target_symbol, t.target_entrez, t.target_ensembl,",
#'                "    i.experiment, i.support_type, i.pubmed_id",
#'                "FROM mirna AS m INNER JOIN mirecords AS i INNER JOIN target",
#'                "AS t ON (m.mature_mirna_uid=i.mature_mirna_uid AND",
#'                "    i.target_uid=t.target_uid)",
#'                "WHERE m.mature_mirna_id='hsa-miR-18a-3p'")
#'   result <- search.multimir(query = qry)
#' 
#' @export search.multimir
search.multimir <- function(query) {
    # To search the multiMiR database on the web server given a MySQL query
    # NOTE: Can only be used after version is set?? due to dbName arg?

    dbName <- getOption("multimir.db.name")
    submit_request(query = query, dbName = dbName, 
                   .cgifields = c("query","dbName"))
}


#' General workhorse function for submitting and returning queries
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get.multimir}.
#' 
#' 
#' @export submit_request
submit_request <- function(url = getOption("multimir.url"), query, ...) {

    request <- RCurl::postForm(url, query = query, ...)
    result  <- XML::readHTMLTable(request)
    parse.multimir(result)

}



#' Parse the Result Returned by the multiMiR Web Server
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get.multimir}.
#' 
#' 
#' @export parse.multimir
parse.multimir <- function(HTML.result) {
    # To parse the result from the multimir web server.  Two tables should
    # return. The first table (result[[1]]) is the summary. And the second table
    # (result[[2]]) has the result in details.

    result <- NULL
    l      <- length(HTML.result)
    if (l == 2) {
        result <- HTML.result[[2]]
    } else if (l == 0) {
        warning(paste("Request to multiMiR web server failed. There could be",
					  "incorrect syntax in your query, or you are not connected",
					  "to the internet.  Alternatively the multiMiR web server",
					  "at http://multimir.ucdenver.edu is temporarily down.\n"))
    } else {
		# This should never happen, but just in case...
		warning(paste("Unexpected result from request of multiMiR web server.",
					  "Problem originates with web server - Please submit issue",
					  "on Github repo"))
	}

    return(result)

}

