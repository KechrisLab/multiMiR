
#' List microRNAs, Genes, Drugs Or Diseases in the multiMiR Package
#' 
#' \code{list.multimir} lists all the unique microRNAs, target genes, drugs, or
#' diseases in the web server of the multiMiR package.
#' 
#' \code{list.multimir} returns all the unique microRNAs, target genes, drugs,
#' or diseases in the web server of the multiMiR package.
#' 
#' @param x a character string indicating what to list. This must be one of the
#' strings \code{"mirna"} (default), \code{"gene"}, \code{"drug"}, or
#' \code{"disease"}. This can be abbreviated and is case insensitive.
#' @param url Deprecated. Use global option \code{multimir.url} instead. 
#' @return \code{list.multimir} returns a data frame with information of
#' microRNAs (microRNA unique ID, organism, mature microRNA accession number,
#' and mature microRNA ID), target genes (gene unique ID, organism, gene
#' symbol, Entrez gene ID, and Ensembl gene ID), drugs (drug names), and
#' diseases (disease name).
#' @author Yuanbin Ru \email{ruyuanbin@@gmail.com}
#' @keywords utilities database
#' @examples
#' 
#'   miRNAs <- list.multimir("mirna")
#'   genes <- list.multimir("gene")
#'   drugs <- list.multimir("drug")
#'   diseases <- list.multimir("disease")
#' 
#' @export list.multimir
list.multimir <- function(x   = c("mirna", "gene", "drug", "disease"),
                          limit = NULL,
                          url = NULL) {

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
    result <- lapply(qry, search.multimir)

    stopifnot(length(result) %in% 1:2)

    # Clean up result and return
    if (length(result) == 2) {
		    stopifnot(names(result[[1]]) == names(result[[2]]))
		    nm      <- names(result[[1]])
        result  <- sort(union(toupper(result[[1]][, 1]), toupper(result[[2]][, 1])))
        result  <- data.frame(result)
		colnames(result) <- nm
    } else {
        result <- result[[1]]
    }

    return(result)

}

