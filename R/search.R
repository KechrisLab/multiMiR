
#' Search the multiMiR Database Given a MySQL Query
#' 
#' This is a function for directly querying the multiMiR database with MySQL
#' queries. Given a MySQL query, it searches and retrieves result from the
#' multiMiR database on the multiMiR web server. To use \code{search_multimir}
#' directly, users will need to be familiar with MySQL and multiMiR table
#' structures. Users are advised to use \code{get_multimir} instead.
#' 
#' search.multimir() has been deprecated and replaced with the search_multimir()
#' version.
#' 
#' @param query a character string for the MySQL query.
#' @return \code{search_multimir} returns a data frame containing results from
#' the multiMiR web server.
#' @keywords utilities database
#' @examples
#' 
#'   ## show all tables in the multiMiR database
#'   tables <- search_multimir(query="show tables")
#' 
#'   ## show the structure of table diana_microt
#'   microt <- search_multimir(query="describe diana_microt")
#' 
#'   ## search for validated target genes of hsa-miR-18a-3p in miRecords
#'   qry <- paste("SELECT m.mature_mirna_acc, m.mature_mirna_id,",
#'                "    t.target_symbol, t.target_entrez, t.target_ensembl,",
#'                "    i.experiment, i.support_type, i.pubmed_id",
#'                "FROM mirna AS m INNER JOIN mirecords AS i INNER JOIN target",
#'                "AS t ON (m.mature_mirna_uid=i.mature_mirna_uid AND",
#'                "    i.target_uid=t.target_uid)",
#'                "WHERE m.mature_mirna_id='hsa-miR-18a-3p'")
#'   result <- search_multimir(query = qry)
#' 
#' @export
search_multimir <- function(query) {
    # To search the multiMiR database on the web server given a MySQL query
    # NOTE: Can only be used after version is set?? due to dbName arg?

    dbName <- getOption("multimir.db.name")
    submit_request(query = query, dbName = dbName, 
                   .cgifields = c("query","dbName"))
}

#' @rdname search_multimir
#' @export
search.multimir <- function(query) {
    .Deprecated("search_multimir")
    search_multimir(query = query)
}




#' General workhorse function for submitting and returning queries
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get_multimir}.
#' 
#' @importFrom XML readHTMLTable
#' @importFrom RCurl postForm
#' @return Table requested in \code{query}.
#' @keywords internal
submit_request <- function(url = full_url("multimir.queries"), query, ...) {

    request <- RCurl::postForm(url, query = query, ... )
    result  <- XML::readHTMLTable(request, stringsAsFactors = FALSE)
    parse_response(result)

}




#' Parse the Result Returned by the multiMiR Web Server
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get_multimir}.
#' 
#' @return The queried table portion of the HTML response.
#' @keywords internal
parse_response <- function(HTML.response) {
    # To parse the response from the multimir web server.  Two tables should
    # return. The first table (response[[1]]) is the summary. And the second
    # table (response[[2]]) has the response in details.

    response <- NULL
    l      <- length(HTML.response)
    if (l == 2) {
        response <- HTML.response[[2]]
    } else if (l > 2) {
        # This should never happen, but just in case...
        stop(paste("Unexpected response from multiMiR web server.",
                   "Problem originates with web server - Please submit issue",
                   "on Github repo"))
    } else if (l == 0) {
        stop("Request to multiMiR web server failed. There could be ",
             "incorrect syntax in your query, or you are not connected to ",
             "the internet.  Alternatively the multiMiR web server at ", 
             "http://multimir.org is temporarily down. \n")
    }  # else if l ==1, just return NULL

    return(response)

}

