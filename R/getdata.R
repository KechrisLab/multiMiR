
#' Get microRNA-target Interactions from the multiMiR Package
#' 
#' The main function to retrieve predicted and validated miRNA-target
#' interactions and their disease and drug associations from the multiMiR
#' package.
#' 
#' get.multimir() has been deprecated and replaced with the get_multimir()
#' version.
#' 
#' \code{get_multimir} is the main and recommended function to retrieve
#' information from the multiMiR package. Input to the function must contain at
#' least one of the followings: miRNA(s), target gene(s), and disease and drug
#' term(s).
#' 
#' The setting of \code{predicted.site} is applicable to three ("miranda",
#' "pita", and "targetscan") of the eight predicted tables.  If
#' \code{predicted.site} is \code{"conserved"}, the function will search
#' conserved target sites annotated by TargetScan, target sites with
#' conservation scores greater than or equal to 0.57 (in human and rat; or
#' 0.566 in mouse) in miRanda, and/or sites with conservation scores greater
#' than or equal to 0.9 in PITA.
#' 
#' Although the summary (if \code{summary=TRUE}) can be used to find results
#' that are recorded by combinations of different databases, please note that
#' for predicted interactions a combination approach may not be as effective as
#' a single algorithm because of age or quality of the tool.
#' 
#' Note: The length of the list supported has been increased from version1.0.1.
#' The size is now limited to 20MB which should accommodate most requests.
#' There is a possibility for technical reasons that the query could fail even
#' if the list is under this limit.  If this occurs it is recommended that you
#' break up the list into smaller batches and submit them sequentially.
#' 
#' @param url Deprecated. The URL for queries is now defined by the package
#' options \code{multimir.url} and \code{multimir.queries}.
#' @param org a character string for the organism. Three organisms are
#' supported so far: human ("hsa" (default), "human", or "Homo Sapiens"), mouse
#' ("mmu", "mouse", or "Mus musculus"), and rat ("rno", "rat", or "Rattus
#' norvegicus"). The organism is case insensitive.
#' @param mirna 'NULL' (default) or a character string or character vector for
#' the mature miRNA(s). It can be the mature miRNA accession number (i.e.
#' "MIMAT0000072"), mature miRNA ID (i.e. "hsa-miR-199a-3p"), or a combination
#' of both (i.e. c("MIMAT0000065", "hsa-miR-30a-5p")).  The character is case
#' insensitive. *See note about the length of list supported.
#' @param target 'NULL' (default) or a character string or character vector for
#' the target gene(s). It can be the gene symbol (i.e. c("TP53", "KRAS")),
#' Entrez gene ID (i.e. c(578, 3845)), Ensembl gene ID (i.e.
#' "ENSG00000171791"), or a combination of any of these identifiers (i.e.
#' c("TP53", 3845, "ENSG00000171791")). The character is case insensitive. *See
#' note about the length of list supported.
#' @param disease.drug 'NULL' (default) or a character string or character
#' vector for the disease(s) and/or drug(s) (i.e. c("bladder cancer",
#' "cisplatin")).  The character is case insensitive.
#' @param table a character string indicating which table(s) in multiMiR to
#' search. Each table contains data from an external database.  Options include
#' "validated" (default, to search all validated tables "mirecords",
#' "mirtarbase", and "tarbase"), "predicted" (to search all predicted tables
#' "diana_microt", "elmmo", "microcosm", "miranda", "mirdb", "pictar", "pita",
#' and "targetscan"), "disease.drug" (to search all disease/drug tables
#' "mir2disease", "pharmaco_mir", and "phenomir"), "all" (to search all of the
#' tables above), or an individual table from above.
#' @param predicted.cutoff.type a character indicating the type of prediction
#' score cutoff. This must be either "p" (default, percentage cutoff) or "n"
#' (number cutoff).
#' @param predicted.cutoff 'NULL' (default) or an integer giving a prediction
#' score cutoff.  By default ('NULL'), the cutoff is '20' (search the top 20\%
#' if \code{predicted.cutoff.type="p"}) or '300000' (search the top 300000 (or
#' all records if total < 300000) if \code{predicted.cutoff.type="n"}).
#' @param predicted.site a character string indicating the type of predicted
#' target sites to search. This can be one of the strings "conserved",
#' "nonconserved", or "all", and can be abbreviated. This only applies to three
#' of the predicted tables ("miranda", "pita", and "targetscan") that have
#' conservation information of the target sites.
#' @param summary logical. Whether to summarize the result (default = FALSE).
#' @param add.link logical. Whether to add link to external database for each
#' result entry.
#' @param use.tibble logical. Whether to use the data_frame class from the
#' tibble package for returned dataframes.  The key benefit for large datasets
#' is more restrictive printing to the console (first 10 rows and only the
#' number of columns that will fit \code{getOption('width')}). See
#' \code{?tible::data_frame} for more information. 
#' @param limit a positive integer. Limits the number of records returned from
#' each table.  Useful in testing potentially large queries.
#' @param legacy.out logical. Whether to return the Bioconductor compatible S4
#' object or the legacy S3 object (default=FALSE).
#'
#' @return \code{get_multimir} returns an S4 object (see
#' \code{?mmquery_bioc-class} containing the queried data and associated
#' metadata. With \code{legacy.out=FALSE} (default), the data is a single
#' dataset with association/interaction type defined by the \code{type}
#' variable. With \code{legacy.out=TRUE} the original S3 object with 3 separate
#' data frames ('predicted', 'validated', and 'disease_drug') is returned. 
#' @keywords utilities database
#' @examples
#'
#'   ## search 'hsa-miR-18a-3p' in validated interactions in human
#'   example1 <- get_multimir(mirna='hsa-miR-18a-3p', summary=TRUE)
#'   columns(example1)
#'   ## target genes that are validated by Luciferase assay
#'   lucif <- select(example1, keytype = "type", keys = "validated", 
#'                   columns = columns(example1))
#'   lucif[grep("Luciferase", lucif$experiment), ]
#'   example1@summary[example1@summary[,"target_symbol"] == "KRAS",]
#'
#'   ## search 'cisplatin' in disease and drug tables in human
#'   example2 <- get_multimir(disease.drug='cisplatin', table='disease.drug')
#'   nrow(example2@data)
#'   head(example2@data)
#'
#' @importFrom purrr map
#' @importFrom tibble as_data_frame
#' @importFrom stats setNames
#' @export get_multimir
get_multimir <- function(url                   = NULL,
                         org                   = "hsa",
                         mirna                 = NULL,
                         target                = NULL,
                         disease.drug          = NULL,
                         table                 = "validated",
                         predicted.cutoff      = NULL,
                         predicted.cutoff.type = "p",
                         predicted.site        = "conserved",
                         summary               = FALSE,
                         add.link              = FALSE,
                         use.tibble            = FALSE,
                         limit                 = NULL,
                         legacy.out            = FALSE) {

    if (!is.null(url)) deprecate_arg("url")
    if (is.null(mirna) & is.null(target) & is.null(disease.drug)) return(NULL) 

    # Argument checking
    if (!table %in% c(all_tables(), "predicted", "validated", "disease.drug",
                      "all")) {
        stop("Invalid table value. See help for options.")
    }
    if (is.null(mirna) & is.null(target) & table == "all") {
        message("Predicted and validated tables require either mirna or ", 
                "target arguments. Only disease/drug tables will be returned.")
        table <- "disease.drug" 
    }

    # Grab argument for storing in return object
    my_args  <- mget(names(formals()), sys.frame(sys.nframe()))
    argmatch <- match(c("table", "org", "mirna", "target", "disease.drug",
                        "predicted.cutoff", "predicted.cutoff.type",
                        "predicted.site"), names(my_args), 0L)
    sqlargs  <- as.list(my_args[c(argmatch)])

    # Parse arguments
    org              <- parse_orgs(org)
    mirna            <- remove_empty_strings(mirna)
    target           <- remove_empty_strings(target)
    predicted.cutoff <- default_cutoff(predicted.cutoff.type, predicted.cutoff)
    .table <- switch(table,
                     all          = all_tables(),
                     validated    = validated_tables(),
                     predicted    = predicted_tables(),
                     disease.drug = diseasedrug_tables(),
                     table)

    # Don't build queries for tables that don't apply to provided arguments
    if (org == "rno") {
        .table <- remove_table(.table, c("diana_microt", "pictar", "pita",
                                         "targetscan"))
    }
    if (is.null(mirna) & is.null(disease.drug)) {
        .table <- remove_table(.table, c("mir2disease", "phenomir"))
    }

    # Build queries and request from server
    queries <- map(.table, build_mmsql, 
                   org                   = org,
                   mirna                 = mirna,
                   target                = target,
                   disease.drug          = disease.drug,
                   predicted.site        = predicted.site,
                   predicted.cutoff.type = predicted.cutoff.type,
                   predicted.cutoff      = predicted.cutoff,
                   limit                 = limit)
    queries   <- setNames(queries, .table)
    # Request data
    .data     <- map(queries, query_multimir, org = org,
                     add.link = add.link, use.tibble = use.tibble)
    # Restructure data and related info for returning
    rtnobject <- extract_mmquery(outlist = .data, org = org, summary = summary,
                                 use.tibble = use.tibble, .args = sqlargs)

    if (add.link) {
        message(paste("Some of the links to external databases may be broken",
                      "due to outdated identifiers in these databases. Please",
                      "refer to Supplementary Table 2 in the multiMiR paper",
                      "for details of the issue.\n"))
    }

    # Choose between S4 and legacy S3 output
    rtnobject <- if (legacy.out) as.mmquery(rtnobject) else as.mmquery_bioc(rtnobject)

    return(rtnobject)

}

#' @rdname get_multimir
#' @export
get.multimir <- function(url                   = NULL,
                         org                   = "hsa",
                         mirna                 = NULL,
                         target                = NULL,
                         disease.drug          = NULL,
                         table                 = "validated",
                         predicted.cutoff      = NULL,
                         predicted.cutoff.type = "p",
                         predicted.site        = "conserved",
                         summary               = FALSE,
                         add.link              = FALSE,
                         use.tibble            = FALSE,
                         limit                 = NULL) {

    .Deprecated("get_multimir")
    get_multimir(url                   = url,
                 org                   = org,
                 mirna                 = mirna,
                 target                = target,
                 disease.drug          = disease.drug,
                 table                 = table,
                 predicted.cutoff      = predicted.cutoff,
                 predicted.cutoff.type = predicted.cutoff.type,
                 predicted.site        = predicted.site,
                 summary               = summary,
                 add.link              = add.link,
                 use.tibble            = use.tibble,
                 limit                 = limit)

}



#' Each org can be specified in one of 3 ways -- this standardizes the argument
#' into the 3 char abbreviation.
#'
#' @return A standardized, abbreviated form of the input org.
#' @keywords internal
parse_orgs <- function(org) {

    # only allows single string (TODO: same as old version?).
    if (!is.null(org)) {
        org <- gsub("hsa|human|homo sapiens", "hsa", org, ignore.case = TRUE)
        org <- gsub("mmu|mouse|mus musculus", "mmu", org, ignore.case = TRUE)
        org <- gsub("rno|rat|rattus norvegicus", "rno", org, ignore.case = TRUE)

        if (!(org %in% c("hsa", "mmu", "rno"))) {
            stop("Organism ", org,  " is not in multiMiR. Current options ",
                 "are 'hsa' (human), 'mmu' (mouse) and 'rno' (rat).\n")
        }     
    } 

    return(org)
}


#' If null, set default predicted.cutoff 
#'
#' @return The default cutoff value. 
#' @keywords internal
default_cutoff <- function(predicted.cutoff.type, predicted.cutoff) {

    if (!is.null(predicted.cutoff.type) & is.null(predicted.cutoff))  {
        predicted.cutoff <- switch(predicted.cutoff.type, 
                                   p = 20, 
                                   n = 300000)
    }

    if (predicted.cutoff.type == "p" & 
        (predicted.cutoff < 1 | predicted.cutoff > 100)) {
        stop(paste("Percent predicted cutoff (predicted.cutoff) should be",
                   "between 1 and 100.\n"))
    }

    return(predicted.cutoff)

}


#' Wrapper for search_multimir for adding feature (printing notification to
#' console)
#' 
#' @return The queried multimir data with the addition of a requested feature.
#' @keywords internal
#' @importFrom tibble as_data_frame
query_multimir <- function(x, org, add.link, use.tibble) {

    cat("Searching", x$table, "...\n")
    x$data <- search_multimir(x$query)

    if (!is.null(x$data)) {
        x$data <- cbind('database' = x$table, x$data, stringsAsFactors = FALSE)
        if (add.link) {
           x$data <- add.multimir.links(x$data, org)
        }
    } else {
        x$data <- data.frame()
    }

    if (use.tibble) x$data <- as_data_frame(x$data)
    if (x$table %in% diseasedrug_tables()) x$data <- unique(x$data)

    return(x)

}


#' Remove tables x from a vector of table names.
#'
#' Typically used when a set of arguments don't apply to a table or would return
#' an error/empty response
#'
#' @param tables A character vector.
#' @param x A second character vector to remove from the first (\code{tables}).
#' @return Character vector \code{tables} excluding the strings matching those
#' in \code{x}.
#' @keywords internal
remove_table <- function(tables, x) tables[!tables %in% x]


#' Remove empty strings from character vector.
#'
#' The WHERE clauses for target and mirna use allow for multiple arguments
#' always separated by 'OR' and several columns are checked for each value 
#' (mirna id, acc; target symbol, entrez, ensemble).  If empty strings "" are
#' present in the get_multimir arguments, Targets and miRNA with empty values in
#' one of these columns will be incorrectly returned.  -- thus purge empty
#' strings first.
#'
#' @param x A character vector
#' @return A character vector with empty strings removed
#' @keywords internal
remove_empty_strings <- function(x) x[x != ""]



