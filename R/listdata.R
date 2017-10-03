
#' List microRNAs, Genes, Drugs Or Diseases in the multiMiR Package
#'
#' \code{list_multimir} lists all the unique microRNAs, target genes, drugs, or
#' diseases in the web server of the multiMiR package.
#'
#' list.multimir() has been deprecated and replaced with the list_multimir()
#' version.
#'
#' @param x a character string indicating what to list. This must be one of the
#' strings \code{"mirna"} (default), \code{"gene"}, \code{"drug"}, or
#' \code{"disease"}. This can be abbreviated and is case insensitive.
#' @param limit a positive integer. Limits the number of records returned from
#' each table.  Useful in testing potentially large queries.
#' @param url Deprecated. Use global option \code{multimir.url} instead.
#' @return \code{list_multimir} returns a data frame with information of
#' microRNAs (microRNA unique ID, organism, mature microRNA accession number,
#' and mature microRNA ID), target genes (gene unique ID, organism, gene
#' symbol, Entrez gene ID, and Ensembl gene ID), drugs (drug names), and
#' diseases (disease name).
#' @author Yuanbin Ru \email{ruyuanbin@@gmail.com}
#' @keywords utilities database
#' @examples
#'   miRNAs <- list_multimir("mirna", limit = 10)
#'   genes <- list_multimir("gene", limit = 10)
#'   drugs <- list_multimir("drug", limit = 10)
#'   diseases <- list_multimir("disease", limit = 10)
#' @importFrom purrr map
#' @export list_multimir
list_multimir <- function(x     = c("mirna", "gene", "drug", "disease"),
                          limit = NULL,
                          url   = NULL) {

    if (!is.null(url)) deprecate_arg("url")
    x   <- match.arg(x)

    # Set chosen query and submit/request from server 
    qry <- switch(x,
                  mirna   = list("SELECT * FROM mirna"),
                  gene    = list("SELECT * FROM target"),
                  drug    = list("SELECT DISTINCT(drug) FROM pharmaco_mir"),
                  disease = list("SELECT DISTINCT(disease) FROM mir2disease",
                                 "SELECT DISTINCT(disease) FROM phenomir"))
    if (!is.null(limit)) qry <- purrr::map(qry, ~ paste(.x, "LIMIT", limit))
    result <- lapply(qry, search_multimir)

    stopifnot(length(result) %in% 1:2)

    # Clean up result and return
    if (length(result) == 2) {
            stopifnot(names(result[[1]]) == names(result[[2]]))
            nm      <- names(result[[1]])
        result  <- sort(union(toupper(result[[1]][, 1]), 
                              toupper(result[[2]][, 1])))
        result  <- data.frame(result)
        colnames(result) <- nm
    } else {
        result <- result[[1]]
    }

    return(result)

}


#' @export list.multimir
#' @rdname list_multimir
list.multimir <- function(x     = c("mirna", "gene", "drug", "disease"),
                          limit = NULL,
                          url   = NULL) {
    .Deprecated("list_multimir")
    list_multimir(x     = x,
                  limit = limit,
                  url   = url)
}


