
#' Get microRNA-target Interactions from the multiMiR Package
#' 
#' The main function to retrieve predicted and validated miRNA-target
#' interactions and their disease and drug associations from the multiMiR
#' package.
#' 
#' \code{get.multimir} is the main and recommended function to retrieve
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
#' 
#' @return \code{get.multimir} returns a list with several data frames
#' containing results from a given external database (e.g., if
#' \code{table="targetscan"}), the predicted (if \code{table= "predicted"}),
#' validated (if \code{table="validated"}), and disease and drug (if
#' \code{table="disease.drug"}) components of multiMiR, and a summary (if
#' \code{summary=TRUE}).
#' @keywords utilities database
#' @examples
#' 
#'   ## search 'hsa-miR-18a-3p' in validated interactions in human
#'   example1 <- get.multimir(mirna='hsa-miR-18a-3p', summary=TRUE)
#'   names(example1)
#'   ## target genes that are validated by Luciferase assay
#'   example1$validated[grep("Luciferase", example1$validated[,"experiment"]),]
#'   example1$summary[example1$summary[,"target_symbol"] == "KRAS",]
#' 
#'   ## search 'cisplatin' in disease and drug tables in human
#'   example2 <- get.multimir(disease.drug='cisplatin', table='disease.drug')
#'   nrow(example2$disease.drug)
#'   head(example2$disease.drug)
#' 
#' @export get.multimir
get.multimir <- function(url = NULL, 
                         org = "hsa", 
                         mirna = NULL, 
                         target = NULL,
                         disease.drug = NULL, 
                         table = "validated",
                         predicted.cutoff = NULL,
                         predicted.cutoff.type = "p",
                         predicted.site = "conserved", 
                         summary = FALSE, 
                         add.link = FALSE) {

    if (!is.null(url)) deprecate_arg("url")
    if (is.null(mirna) & is.null(target) & is.null(disease.drug)) return(NULL) 

    # Don't use scientific notation in fn environment 
    #   for converting num to char where "3e4" != "30000"
    scipen.orig <- getOption("scipen")
    options(scipen = 999)
    on.exit(options(scipen = scipen.orig))

    # Collect args for use in query builders
    my_args  <- mget(names(formals()), sys.frame(sys.nframe()))
    argmatch <- match(c("org", "mirna", "target", "disease.drug",
                        "predicted.cutoff", "predicted.cutoff.type",
                        "predicted.site"), names(my_args), 0L)
    sqlargs  <- as.list(my_args[c(argmatch)])

    # Argument parsing and checking
    check_tables(table)
    org              <- parse_orgs(org)
    predicted.cutoff <- default_cutoff(predicted.cutoff.type, predicted.cutoff)

    # Other default table names and reassign changed args
    sqlargs$mirna.table      <- "mirna"
    sqlargs$target.table     <- "target"
    sqlargs$org              <- org
    sqlargs$mirna            <- mirna        
    sqlargs$target           <- target       
    sqlargs$disease.drug     <- disease.drug 
    sqlargs$predicted.cutoff <- predicted.cutoff

#     get_query <- function(x, ...) {
#         this_query <- do.call(x[["query_name"]], c(table = x[["table"]], ...))
#         c(x, "query" = this_query)
#         
#     }

    # NOTE: Testing version
    get_query <- function(x, ...) {
        this_query <- do.call(x[["query_name"]], c(table = x[["table"]], ...))
        #c(x, "query" = this_query)
        this_query
    }

    tbls_to_query <- table_query_lookup(table, mirna, target)
    my_queries    <- lapply(tbls_to_query, get_query, sqlargs)
    names(my_queries) <- do.call(c, lapply(tbls_to_query, function(x) x$table))

    # NOTE: Commented out for testing -- comparing queries w/ old vers
#     result <- lapply(my_queries, function(x) {
#                          cat("Searching", x$table, "...\n")
#                          search.multimir(x$query)
#                                            })
#     return(list(.query = my_queries, .data = result))
    # NOTE: Restructuring for testing

    return(my_queries)


    # TODO:
    # 1) add table name to each dataset (with varname='database'), 
    # 2) rbind all requested tables, 
    # 3) if add.link add link  
    # 4) if summary add summary.
    # 5) set type, query_name, and table as attributes

#    if (add.link & !is.null(result[[table]])) 
#        result[[table]] <- add.multimir.links(result[[table]], org)
#    }
#    if (summary) {
#        result[["summary"]] <- multimir.summary(result)
#    }
#    return(result)
}



#' Create/filter lookup table for all SQL tables to query
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get.multimir}.
#'
#' @keywords internal
table_query_lookup <- function(tbl_arg, mirna, target) {

    factor_op <- getOption("stringsAsFactors")
    options(stringsAsFactors = FALSE)
    on.exit(options(stringsAsFactors = factor_op))

    table_query_lookup <- 
        rbind(data.frame(type       = c("validated"), 
                         query_name = c("query_validated"), 
                         table      = c("mirecords", "mirtarbase", "tarbase")),
              data.frame(type       = c("predicted"), 
                         query_name = c("query_predicted"), 
                         table      = c("diana_microt", "elmmo", "microcosm",
                                        "miranda", "mirdb", "pictar", "pita",
                                        "targetscan")),
              data.frame(type       = c("disease.drug"), 
                         query_name = c("query_disease"),
                         table      = c("mir2disease", "pharmaco_mir",
                                        "phenomir")))

    tbl_arg <- tolower(tbl_arg)
    # Either mirna or target required for predicted and validated tables 
    if (is.null(mirna) & is.null(target) & tbl_arg == "all") {
        message("Predicted and validated tables require either mirna or target ",
                "arguments. Only disease/drug tables will be returned.")
        tbl_arg <- "disease.drug" 
    }

	# Create list of tables to query (input value has to be length 1L, handled above)
	tables <- switch(tbl_arg,
                     validated    = c("mirecords", "mirtarbase", "tarbase"),
                     predicted    = c("diana_microt", "elmmo", "microcosm",
                                      "miranda", "mirdb", "pictar", "pita",
                                      "targetscan"),
                     disease.drug = c("mir2disease", "pharmaco_mir", "phenomir"),
                     all 		  = c("mirecords", "mirtarbase", "tarbase",
                                      "diana_microt", "elmmo", "microcosm",
                                      "miranda", "mirdb", "pictar", "pita",
                                      "targetscan", "mir2disease", "pharmaco_mir",
                                      "phenomir"),
                     tolower(tbl_arg))

    rtn <- subset(table_query_lookup, table %in% tables)
    rtn <- split(rtn, seq(nrow(rtn)))
    return(rtn)

}


#' Each org can be specified in one of 3 ways -- this standardizes the argument
#' into the 3 char abbreviation.
#'
#' @keywords internal
parse_orgs <- function(org) {

	# Currently only allows single string (TODO: same as old version?).
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




#' Check table argument matches valid table name
#'
#' @keywords internal
check_tables <- function(table, force_update = FALSE) {

	stopifnot(!is.null(table), length(table) == 1L) 

    tbl_lst <- getOption("multimir.tables.list")

    if (force_update | is.null(tbl_lst)) tbl_lst <- multimir_dbTables()

    valid_tbls <- c(setdiff(tbl_lst, c("map_counts", "map_metadata", "metadata")), 
                    "validated", "predicted", "disease.drug", "all")

    if (!table %in% valid_tbls) {
        stop(paste("Table", table, "does not exist!\n", "Please use",
                   "'multimir_dbTables()' to see a list of available",
                   "tables.\n"))
    } else {
        invisible()
    }

}


#' If null, set default predicted.cutoff 
#'
#' @keywords internal
default_cutoff <- function(predicted.cutoff.type, predicted.cutoff) {

    if (!is.null(predicted.cutoff.type) & is.null(predicted.cutoff))  {
        predicted.cutoff <- switch(predicted.cutoff.type, 
                                   p = 20, 
                                   n = 300000)
	}

	if (predicted.cutoff.type == "p" & (predicted.cutoff < 1 | predicted.cutoff > 100)) {
		stop(paste("Percent predicted cutoff (predicted.cutoff) should be",
				   "between 1 and 100.\n"))
	}

    return(predicted.cutoff)

}


