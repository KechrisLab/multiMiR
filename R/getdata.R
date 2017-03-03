
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
#' @return \code{get.multimir} returns a list with several data frames
#' containing results from a given external database (e.g., if
#' \code{table="targetscan"}), the predicted (if \code{table= "predicted"}),
#' validated (if \code{table="validated"}), and disease and drug (if
#' \code{table="disease.drug"}) components of multiMiR, and a summary (if
#' \code{summary=TRUE}).
#' @author Yuanbin Ru \email{ruyuanbin@@gmail.com}
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

    # Don't use scientific notation in fn environment 
    #   for converting num to char where "3e4" != "30000"
    scipen.orig <- getOption("scipen")
    options(scipen = 999)
    on.exit( options(scipen = scipen.orig) )

    # Collect args for use in query builders
    my_args  <- mget(names(formals()), sys.frame(sys.nframe()))
    argmatch <- match(c("org", "mirna", "target", "disease.drug",
                        "predicted.cutoff", "predicted.cutoff.type",
                        "predicted.site"), names(my_args), 0L)
    sqlargs  <- as.list(my_args[c(argmatch)])

    # Other default table names
    sqlargs$mirna.table  <- "mirna"
    sqlargs$target.table <- "target"

    # Prep table argument for inclusion in sql query
	stopifnot(!is.null(table), length(table) == 1L) 
    multimir.tables <- c(setdiff(as.character(multimir_dbTables()),
                                 c("map_counts", "map_metadata", "metadata")),
                         "validated", "predicted", "disease.drug", "all")
    if (!table %in% multimir.tables) {
        stop(paste("Table", table, "does not exist!\n", "Please use",
                   "'multimir_dbTables()' to see a list of available",
                   "tables.\n"))
    } 

    if (is.null(mirna) & is.null(target) & is.null(disease.drug)) return(NULL) 

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

    # Prep these names for use in SQL query
    wrap_in_parens <- function(x) {
        if (!is.null(x)) paste0("('", paste(x, collapse = "','"), "')")
    }
    mirna        <- wrap_in_parens(mirna)
    target       <- wrap_in_parens(target)
    disease.drug <- wrap_in_parens(disease.drug)

	# Set default predicted.cutoff 
    if (!is.null(predicted.cutoff.type) & is.null(predicted.cutoff))  {
        predicted.cutoff <- switch(predicted.cutoff.type, 
                                           p = 20, n = 300000)
	}
	if (predicted.cutoff.type == "p" & (predicted.cutoff < 1 | predicted.cutoff > 100)) {
		stop(paste("Percent predicted cutoff (predicted.cutoff) should be",
				   "between 1 and 100.\n"))
	}

    # Reassign changed args
    sqlargs$org              <- org
    sqlargs$mirna            <- mirna        
    sqlargs$target           <- target       
    sqlargs$disease.drug     <- disease.drug 
    sqlargs$predicted.cutoff <- predicted.cutoff

    get_query <- function(x, ...) {
        cl <- do.call(x[["query_name"]], c(table = x[["table"]]), a.list(...))
    }

    cat("Searching", table, "...\n")
    tbls_to_query <- table_query_lookup(table)
    my_queries    <- apply(tbls_to_query, 1, get_query, sqlargs)

    return(my_queries)

#	rtn <- lapply(c("mirecords", "mirtarbase", "tarbase"), get.multimir.by.table, 
#				  url = url, org = org, mirna = mirna, target = target)
#	do.call(rbind, rtn)
#
#	if (add.link & !is.null(result[[table]])) {
#            result[[table]] <- add.multimir.links(result[[table]], 
#                                                  org)
#        }
#
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
#' @param tbl_arg PLACEHOLDER
#' 
table_query_lookup <- function(tbl_arg) {

    factor_op <- getOption("stringsAsFactors")
    options(stringsAsFactors = FALSE)
    on.exit(options(stringsAsFactors = factor_op))

	# Create list of tables to query (input value has to be length 1L, handled above)
    table_query_lookup <- 
        rbind(data.frame(type = "validated", query_name = "query_validated", 
                         table = c("mirecords", "mirtarbase", "tarbase")),
              data.frame(type = "predicted", query_name = "query_predicted", 
                         table = c("diana_microt", "elmmo", "microcosm",
                                   "miranda", "mirdb", "pictar", "pita",
                                   "targetscan")),
              data.frame(type = "disease.drug", query_name = "query_disease",
                         table = c("mir2disease", "pharmaco_mir", "phenomir")))

	tables <- 
		switch(tolower(tbl_arg),
			   validated    = c("mirecords", "mirtarbase", "tarbase"),
			   predicted    = c("diana_microt", "elmmo", "microcosm",
								"miranda", "mirdb", "pictar", "pita",
								"targetscan"),
			   disease.drug = c("mir2disease", "pharmaco_mir", "phenomir"),
			   all 		    = c("mirecords", "mirtarbase", "tarbase",
								"diana_microt", "elmmo", "microcosm",
								"miranda", "mirdb", "pictar", "pita",
								"targetscan", "mir2disease", "pharmaco_mir",
								"phenomir"),
			   tolower(tbl_arg))

    rtn <- subset(table_query_lookup, table %in% tables)
    return(rtn)

}
 
 
 
# #' Get microRNA/target Information from a Given Table in the multiMiR Package
# #' 
# #' This is an internal multiMiR function that is not intended to be used
# #' directly.  Please use \code{get.multimir}.
# #'
# # @param url                   PLACEHOLDER
# # @param org                   PLACEHOLDER
# # @param mirna                 PLACEHOLDER
# # @param target                PLACEHOLDER
# # @param table                 PLACEHOLDER
# # @param disease.drug          PLACEHOLDER
# # @param predicted.cutoff      PLACEHOLDER
# # @param predicted.cutoff.type PLACEHOLDER
# # @param predicted.site        PLACEHOLDER
# # @param mirna.table           PLACEHOLDER
# # @param target.table          PLACEHOLDER
# get.multimir.by.table <- function(table                 = NULL,
#                                   url                   = NULL,
#                                   org                   = c("hsa", "mmu", "rno"),
#                                   mirna                 = NULL,
#                                   target                = NULL,
#                                   disease.drug          = NULL,
#                                   predicted.cutoff      = NULL,
#                                   predicted.cutoff.type = "p",
#                                   predicted.site        = "conserved",
#                                   mirna.table           = "mirna",
#                                   target.table          = "target") {
#     # To get miRNA-target interactions in a given table
#     
#     # prepare query for validated target table
# 	if (table %in% c("mirecords", "mirtarbase", "tarbase")) {
# 		qry <- query_validated(table        = table,
#                                mirna        = mirna,
#                                target       = target,
# 							   org          = org,
#                                mirna.table  = mirna.table,
# 							   target.table = target.table)
# 	}
# 
# 	get.multimir.by.table(table                 = table,
#                           org                   = org,
# 						  mirna                 = mirna,
# 						  target                = target,
# 						  disease.drug          = disease.drug,
# 						  predicted.cutoff      = predicted.cutoff,
# 						  predicted.cutoff.type = predicted.cutoff.type,
# 						  predicted.site        = predicted.site)
#     
#     # prepare query for predicted target table
# 	if (table %in% c("diana_microt", "elmmo", "microcosm", "miranda", "mirdb",
# 					 "pictar", "pita", "targetscan")) {
# 		qry <- query_predicted(table                 = table,
#                                mirna                 = mirna,
#                                target                = target,
# 							   org                   = org,
#                                mirna.table           = mirna.table,
# 							   target.table          = target.table,
#                                predicted.cutoff      = predicted.cutoff,
#                                predicted.cutoff.type = predicted.cutoff.type,
#                                predicted.site        = predicted.site)
# 	}
# 
#     
#     # prepare query for miRNA-disease/drug table
# 	if (table %in% c("mir2disease", "pharmaco_mir", "phenomir")) {
# 
# 		qry <- query_disease(table        = table,
#                              mirna        = mirna,
#                              target       = target,
#                              org          = org,
# 							 mirna.table  = mirna.table,
#                              target.table = target.table,
# 							 disease.drug = disease.drug)
# 	}
# 
#     
#     # query the database
#     result <- search.multimir(query = qry)
#     if (!is.null(result)) 
#         result <- cbind(database = table, result)
#     if (table %in% c("mir2disease", "pharmaco_mir", "phenomir")) 
#         result <- unique(result)
# 
#     return(result)
# }
# 
# 

#' Prepare query for validated target table
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get.multimir}.
#'
#' @param table One of the three validated tables
#' @keywords internal
query_validated <- function(table = c("mirecords", "mirtarbase", "tarbase"), 
                            mirna = NULL, target = NULL, org = NULL,
                            mirna.table, target.table, ...) {

    if (is.null(mirna) & is.null(target)) {
        return(NULL)
    }
	table <- match.arg(table)
    paste("SELECT m.mature_mirna_acc, m.mature_mirna_id,",
          "t.target_symbol, t.target_entrez, t.target_ensembl,",
          "i.experiment, i.support_type, i.pubmed_id",
          "FROM", mirna.table, "AS m INNER JOIN", table, 
          "AS i INNER JOIN", target.table, 
          "AS t ON (m.mature_mirna_uid=i.mature_mirna_uid AND",
          "i.target_uid=t.target_uid) WHERE",
          subquery_mirnatarget(mirna, target), 
          "AND",
          subquery_org(table, org))

}


#' Prepare query for predicted target table
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get.multimir}.
#'
query_predicted <- function(table = c("diana_microt", "elmmo", "microcosm",
                                      "miranda", "mirdb", "pictar", "pita",
                                      "targetscan"), 
                            mirna, target, org, mirna.table, target.table,
                            predicted.cutoff, predicted.cutoff.type,
                            predicted.site = c("conserved", "nonconserved", "all"),
                            ...) {
	table <- match.arg(table)
	predicted.site <- match.arg(predicted.site)

    if (is.null(mirna) & is.null(target)) return(NULL)
    # NOTE: Should this be a message saying org 'rno' doesn't exist for these tables?
    if (org == "rno" & table %in% c("diana_microt", "pictar", "pita",
                                    "targetscan")) return(NULL)

    # Set name of score variable by table
    score_vars <- switch(table,
                         diana_microt = "i.miTG_score",
                         elmmo        = "i.p ",
                         microcosm    = "i.score", # NOTE: didn't have the 'AS score' text
                         mirdb        = "i.score", #       didn't have the 'AS score' text
                         pictar       = "i.score", #       didn't have the 'AS score' text
                         miranda      = "i.mirsvr_score",
                         pita         = "i.ddG",
                         targetscan   = "i.context_plus_score")

    # Create qry_conserve and 'name' variable 
    conserved_list <- subquery_conserved(predicted.site = predicted.site, 
                                        table = table, org = org)


    # Basic SQL structure:
    # predicted: X(choose score var) + Xbase + X(mirna and/or target) + Xorg + #
    #            Xconserved(for 3 of 8) + Xcutoff
	qry_list <- 
        list(base = paste("SELECT m.mature_mirna_acc, m.mature_mirna_id,",
                          "t.target_symbol, t.target_entrez, t.target_ensembl,",
                          score_vars,
                          "AS score FROM", mirna.table, "AS m INNER JOIN", table, 
                          "AS i INNER JOIN", target.table, 
                          "AS t ON (m.mature_mirna_uid=i.mature_mirna_uid",
                          "AND i.target_uid=t.target_uid) WHERE"),
             mirna_target = subquery_mirnatarget(mirna, target),
             org          = subquery_org(table = table, org = org),
             conserved    = conserved_list$query,
             cutoff       = subquery_cutoff(table, predicted.cutoff.type, 
                                            predicted.cutoff,
                                            conserve_list$name, score_vars))

    qry_list <- qry_list[do.call(c, lapply(qry_list, function(x) !is.na(x)))]
    qry      <- paste(qry_list[[1]], paste(do.call(c, qry_list[-1]), collapse = " AND ")) 

	return(qry)
}


#' Prepare query for miRNA-disease/drug table
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get.multimir}.
#'
query_disease <- function(table = c("mir2disease", "pharmaco_mir", "phenomir"), 
						  mirna, target, org, mirna.table, target.table,
						  disease.drug, ...) {

	table <- match.arg(table)
    if (table %in% c("mir2disease", "phenomir") & 
        is.null(mirna) & is.null(disease.drug)) return(NULL)

    qry_base <- 
        switch(table,
               mir2disease = paste("SELECT m.mature_mirna_acc, m.mature_mirna_id, 'NA' AS",
                                   "target_symbol, 'NA' AS target_entrez, 'NA' AS",
                                   "target_ensembl, i.disease AS disease_drug,",
                                   "CONCAT_WS('. ', i.year, i.title)",
                                   "AS paper_pubmedID FROM", 
                                   mirna.table, "AS m INNER JOIN", table, 
                                   "AS i ON (m.mature_mirna_uid=i.mature_mirna_uid) WHERE"),
               pharmaco_mir = paste("SELECT m.mature_mirna_acc, m.mature_mirna_id,",
                                    "t.target_symbol, t.target_entrez, t.target_ensembl,",
                                    "i.drug AS disease_drug, i.pubmed_id AS paper_pubmedID",
                                    "FROM", mirna.table, "AS m INNER JOIN", table, 
                                    "AS i INNER JOIN", target.table, 
                                    "AS t ON (m.mature_mirna_uid=i.mature_mirna_uid AND",
                                    "i.target_uid=t.target_uid) WHERE"),
               phenomir     = paste("SELECT m.mature_mirna_acc, m.mature_mirna_id, 'NA' AS",
                                    "target_symbol, 'NA' AS target_entrez, 'NA' AS",
                                    "target_ensembl, i.disease AS disease_drug, i.pubmed_id AS",
                                    "paper_pubmedID FROM", mirna.table, "AS m INNER JOIN", table, 
                                    "AS i ON (m.mature_mirna_uid=i.mature_mirna_uid) WHERE"))
    qry_diseasedrug <- switch(table, 
                              mir2disease  = paste("i.disease IN", disease.drug),
                              pharmaco_mir = paste("i.drug IN", disease.drug),
                              phenomir     = paste("(i.disease IN", disease.drug, "OR",
                                                   "i.disease_class IN", disease.drug, ")"))
    qry_diseasedrug <- ifelse(is.null(disease.drug), "", qry_diseasedrug)

    qry_org    <- subquery_org(table = table, org = org)
    qry_target <- ifelse(table != "pharmaco_mir", "",
                         subquery_mirnatarget(target = target))
    qry_mirna  <- subquery_mirnatarget(mirna = mirna)
    
    # Basic SQL structure is:
    #   mir2disease: base + disease.drug + org + mirna
    #   pharmaco_mir: base + target + disease.drug + org + mirna
    #   phenomir: base +  disease.drug + org + mirna
    qry_list <- list(qry_base, qry_target, qry_diseasedrug, qry_org, qry_mirna)
    qry_list <- qry_list[do.call(c, lapply(qry_list, function(x) x != ""))]
    qry      <- paste(qry_list[[1]], paste(do.call(c, qry_list[-1]), collapse = " AND ")) 

	return(qry)
}



#' Internal function for building mirna/target portion of query
#'
#' @param mirna see ?get.multimir
#' @param target see ?get.multimir
#' @keywords internal
subquery_mirnatarget <- function(mirna = NULL, target = NULL) {

	stopifnot(!(is.null(mirna) & is.null(target)))

    qry_list <- 
        list(mirna  = paste("( m.mature_mirna_acc IN", mirna, 
                            "OR m.mature_mirna_id IN", mirna, ")"),
             target = paste("( t.target_symbol IN", target, 
                            "OR t.target_entrez IN", target, 
                            "OR t.target_ensembl IN", target, ")"))
    qry_logic <- 
        list(mirna          = !is.null(mirna),
             target         = !is.null(target))

	qry_list <- qry_list[do.call(c, qry_logic)]
	do.call(paste, c(qry_list, sep = " AND "))

}

#' Internal function for building org portion of query
#' 
#' @param table One of the validate, predict, or disease tables
#' @param org One of the 3 organisms (hsa, mmu, rno)
#' @keywords internal
subquery_org <- function(table, org) {

    stopifnot(!is.null(org), !is.null(table))

    if (table %in% c("tarbase", "mir2disease", "phenomir")) {
        paste("m.org =", quote_wrap(org))
    } else {
		paste("m.org =", quote_wrap(org), "AND t.org =", quote_wrap(org))
    }

}

#' Internal function for building conserved (predicted.site) portion of query
#' 
#' Process conserved/nonconserved sites in miranda, pita & targetscan
#' (other tables don't have conserved/nonconserved annotation)
#'   - 2 objects created: 'name' and qry_conserve
#'   - qry_conserv structure: conserved_var + operator + cut_value
#'
#' @param predicted.site see ?get.multimir
#' @param table see ?get.multimir
#' @param org see ?get.multimir
#' @keywords internal
subquery_conserved <- function(predicted.site, table, org) {

    if (table %in% c("miranda", "pita", "targetscan")) {
        conserv_args <- 
            list(suffix    = switch(predicted.site,
                                    conserved    = "c1",
                                    nonconserved = "c0"),
                 variable  = switch(table,
                                    pita = "i.conservation",
                                    miranda = "i.conservation",
                                    targetscan = "i.conserved_site"),
                 operator  = switch(table,
                                    pita       = ifelse(predicted.site == "conserved", ">=", "<"),
                                    miranda    = ifelse(predicted.site == "conserved", ">=", "<"),
                                    targetscan = "="),
                 cut_value = switch(table,
                                    pita       = "0.9",
                                    miranda    = ifelse(org == "mmu", "0.566", "0.57"),
                                    targetscan = "'N'"))

        name <- paste(table, org, conserv_args$suffix, sep = ".")
        qry  <- paste(do.call(c, conserv_args[-1]), collapse = " ")
    } else {
        name <- paste(table, org, sep = ".")
        qry  <- ""
    }

    return(list(`name` = name, `query` = qry))

}

#' Internal function for creating the portion of the query for filtering based
#' on score
#' 
#' Depends on query_predicted() for input arguments.
#'
#' @param table
#' @param predicted.cutoff.type
#' @param predicted.cutoff
#' @param name
#' @param score_vars
#' @keywords internal
subquery_cutoff <- function(table, predicted.cutoff.type, predicted.cutoff,
                            name, score_vars) {

	cutoffs          <- get.multimir.cutoffs()
    tbl_count        <- cutoffs[[name]][["count"]]
    count_min        <- 10000
    count_max        <- 300000
    cutoff_too_small <- paste("Number predicted cutoff (predicted.cutoff)",
                              predicted.cutoff, "may be too small. A cutoff",
                              "of 10000 will be used instead.\n")
    cutoff_too_large <- paste0("Number predicted cutoff (predicted.cutoff) ",
                               predicted.cutoff, " is larger than the total ",
                               "number of records in table ", table, 
                               ". All records will be queried.\n")

	# get dataset-specific score cutoff
	if (predicted.cutoff.type == "p") {
        score.cutoff <- cutoffs[[name]][[paste0(predicted.cutoff, "%")]]
	} else if (predicted.cutoff.type == "n") {
        if (predicted.cutoff < count_min) message(cutoff_too_small)
        if (predicted.cutoff > tbl_count) message(cutoff_too_large)

        cut_pred.cutoff <- cut(predicted.cutoff, labels = FALSE,
                               breaks = c(0, count_min, tbl_count))
        adj.pred.cutoff <- switch(as.character(cut_pred.cutoff), 
                                  '1' = count_min, '2' = predicted.cutoff, 
                                  '3' = tbl_count, NA)
        adj.pred.cutoff <- ifelse(tbl_count < adj.pred.cutoff & 
                                  tbl_count >= count_min, 
                                  tbl_count, adj.pred.cutoff)
        score.cutoff    <- ifelse(is.na(adj.pred.cutoff), NA,
                                  cutoffs[[name]][[as.character(adj.pred.cutoff)]])
    }

    operator <- ifelse(table %in% c("miranda", "pita", "targetscan"), "<=", ">=")
    qry      <- ifelse(is.na(score.cutoff), NA, 
                       paste(score_vars, ">=", score.cutoff, "ORDER BY",
                             score_vars, "DESC"))
	return(qry)

}


# #' Get Validated microRNA-target Interactions from the multiMiR Package
# #' 
# #' This is an internal multiMiR function that is not intended to be used
# #' directly.  Please use \code{get.multimir}.
# #' 
# # @examples
# # 	 get.multimir.validated(mirna = "hsa-miR-18a-3p")
# # 
# get.multimir.validated <- function(url = NULL, org = "hsa", 
#                                    mirna = NULL, target = NULL) {
# 
#     if (!is.null(url)) deprecate_arg("url")
#     if (is.null(mirna) & is.null(target)) return(NULL)
# 
#     # get.multimir.by.table for each table in c() 
# 	rtn <- lapply(c("mirecords", "mirtarbase", "tarbase"), get.multimir.by.table, 
# 				  url = url, org = org, mirna = mirna, target = target)
# 	do.call(rbind, rtn)
# 
# }
# 
# #' Get Predicted microRNA-target Interactions from the multiMiR Package
# #' 
# #' This is an internal multiMiR function that is not intended to be used
# #' directly.  Please use \code{get.multimir}.
# #' 
# get.multimir.predicted <- function(url = NULL, org = "hsa", mirna = NULL, 
#                                    target = NULL, cutoff = NULL, 
# 								   cutoff.type = "p", site = "conserved") {
# 
#     if (!is.null(url)) deprecate_arg("url")
#     if (is.null(mirna) & is.null(target)) return(NULL)
# 
#     # get.multimir.by.table for each table in c() 
# 	rtn <- lapply(c("diana_microt", "elmmo", "microcosm", "miranda", "mirdb",
# 					"pictar", "pita", "targetscan"),
# 				  get.multimir.by.table, org = org, mirna = mirna, 
# 				  target = target, predicted.cutoff = cutoff,
# 				  predicted.cutoff.type = cutoff.type, predicted.site = site)
# 	do.call(rbind, rtn)
# 
# }
# 
# #' Get microRNA-disease/drug interactions from the multiMiR Package
# #' 
# #' This is an internal multiMiR function that is not intended to be used
# #' directly.  Please use \code{get.multimir}.
# #' 
# get.multimir.disease <- function(url = NULL,
#                                  org = "hsa",
#                                  mirna = NULL,
#                                  target = NULL,
#                                  disease.drug = NULL) {
# 
#     if (!is.null(url)) deprecate_arg("url")
# 	if (is.null(mirna) & is.null(target) & is.null(disease.drug)) return(NULL)
# 
#     # get.multimir.by.table for each table in c() 
# 	rtn <- lapply(c("mir2disease", "pharmaco_mir", "phenomir"),
# 				  get.multimir.by.table, org = org, mirna = mirna, 
# 				  target = target, disease.drug = disease.drug)
# 	do.call(rbind, rtn)
# 
# }
