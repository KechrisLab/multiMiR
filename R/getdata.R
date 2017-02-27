
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
#' @param url a character string for the URL of the multiMiR web server.  The
#' default is getOption(" multimir.url")
#' ("http://multimir.ucdenver.edu/cgi-bin/multimir.pl").
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
get.multimir <- function(url = getOption("multimir.url"), 
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
    # The main function to search miRNA-target and miRNA-disease interactions
    result <- list()
    if (table %in% c("all", "validated", "predicted", "disease.drug")) {
        if (table == "predicted" | table == "all") {
            # search predicted miRNA-target tables
            result[["predicted"]] <- 
                get.multimir.predicted(url = url, 
                                       org = org, 
                                       mirna = mirna,
                                       target = target, 
                                       cutoff = predicted.cutoff, 
                                       cutoff.type = predicted.cutoff.type, 
                                       site = predicted.site)
            if (add.link & !is.null(result[["predicted"]])) {
                result[["predicted"]] <- 
                    add.multimir.links(result[["predicted"]], org)
            }
        }
        if (table == "validated" | table == "all") {
            # search validated miRNA-target tables
            result[["validated"]] <- 
                get.multimir.validated(url = url, 
                                       org = org, 
                                       mirna = mirna, 
                                       target = target)
            if (add.link & !is.null(result[["validated"]])) {
                result[["validated"]] <- add.multimir.links(result[["validated"]], 
                                                            org)
            }
        }
        if (table == "disease.drug" | table == "all") {
            # search miRNA-disease tables
            result[["disease.drug"]] <- 
                get.multimir.disease(url = url,
                                     org = org,
                                     mirna = mirna,
                                     target = target,
                                     disease.drug = disease.drug)
            if (add.link & !is.null(result[["disease.drug"]])) {
                result[["disease.drug"]] <- add.multimir.links(result[["disease.drug"]], 
                                                               org)
            }
        }
    } else {
        # search an individual table
        result[[table]] <- get.multimir.by.table(url = url,
                                                 org = org,
                                                 mirna = mirna,
                                                 target = target,
                                                 disease.drug = disease.drug,
                                                 table = table,
                                                 predicted.cutoff = predicted.cutoff,
                                                 predicted.cutoff.type = predicted.cutoff.type,
                                                 predicted.site = predicted.site)
        if (add.link & !is.null(result[[table]])) {
            result[[table]] <- add.multimir.links(result[[table]], 
                                                  org)
        }
    }

    if (summary) {
        result[["summary"]] <- multimir.summary(result)
    }
    return(result)
}


# To get validated miRNA-target interactions


#' Get Validated microRNA-target Interactions from the multiMiR Package
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get.multimir}.
#' 
#' 
#' @export get.multimir.validated
get.multimir.validated <- function(url = NULL,
                                   org = "hsa", 
                                   mirna = NULL, 
                                   target = NULL) {
    if (is.null(mirna) & is.null(target)) {
        return(NULL)
    }
    result <- NULL
    for (table in c("mirecords", "mirtarbase", "tarbase")) {
        r <- get.multimir.by.table(url = url, 
                                   org = org, 
                                   mirna = mirna, 
                                   target = target, 
                                   table = table)
        result <- rbind(result, r)
    }
    return(result)
}


# To get predicted miRNA-target interactions


#' Get Predicted microRNA-target Interactions from the multiMiR Package
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get.multimir}.
#' 
#' 
#' @export get.multimir.predicted
get.multimir.predicted <- function(url = NULL, 
                                   org = "hsa", 
                                   mirna = NULL, 
                                   target = NULL, 
                                   cutoff = NULL, 
                                   cutoff.type = "p", 
                                   site = "conserved") {
    if (is.null(mirna) & is.null(target)) {
        return(NULL)
    }
    result <- NULL
    for (table in c("diana_microt", "elmmo", "microcosm", "miranda", "mirdb",
                    "pictar", "pita", "targetscan")) {
        r <- get.multimir.by.table(url = url, 
                                   org = org, 
                                   mirna = mirna, 
                                   target = target, 
                                   table = table, 
                                   predicted.cutoff = cutoff, 
                                   predicted.cutoff.type = cutoff.type, 
                                   predicted.site = site)
        result <- rbind(result, r)
    }
    return(result)
}


# To get miRNA-disease/drug interactions


#' Get microRNA-disease/drug Relationships from the multiMiR Package
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get.multimir}.
#' 
#' 
#' @export get.multimir.disease
get.multimir.disease <- function(url = NULL,
                                 org = "hsa",
                                 mirna = NULL,
                                 target = NULL,
                                 disease.drug = NULL) {
    if (is.null(mirna) & is.null(target) & is.null(disease.drug)) {
        return(NULL)
    }
    result <- NULL
    for (table in c("mir2disease", "pharmaco_mir", "phenomir")) {
        r <- get.multimir.by.table(url = url,
                                   org = org,
                                   mirna = mirna,
                                   target = target,
                                   table = table,
                                   disease.drug = disease.drug)
        result <- rbind(result, r)
    }
    return(result)
}



#' Get microRNA/target Information from a Given Table in the multiMiR Package
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get.multimir}.
#'
#' @param url                   PLACEHOLDER
#' @param org                   PLACEHOLDER
#' @param mirna                 PLACEHOLDER
#' @param target                PLACEHOLDER
#' @param table                 PLACEHOLDER
#' @param disease.drug          PLACEHOLDER
#' @param predicted.cutoff      PLACEHOLDER
#' @param predicted.cutoff.type PLACEHOLDER
#' @param predicted.site        PLACEHOLDER
#' @param mirna.table           PLACEHOLDER
#' @param target.table          PLACEHOLDER
#' @export get.multimir.by.table
get.multimir.by.table <- function(url                   = NULL,
                                  org                   = "hsa",
                                  mirna                 = NULL,
                                  target                = NULL,
                                  table                 = NULL,
                                  disease.drug          = NULL,
                                  predicted.cutoff      = NULL,
                                  predicted.cutoff.type = "p",
                                  predicted.site        = "conserved",
                                  mirna.table           = "mirna",
                                  target.table          = "target") {
    # To get miRNA-target interactions in a given table

    multimir.tables <- setdiff(as.character(multimir_dbTables()),
                               c("map_counts", "map_metadata", "metadata"))
    if (!table %in% multimir.tables) {
        stop(paste("Table", table, "does not exist!\n", "Please use",
                   "'multimir_dbTables()' to see a list of available tables.\n", 
            sep = " "))
    }
    cat("Searching", table, "...\n")
    
    if ((is.null(mirna) & is.null(target) & is.null(disease.drug)) |
        is.null(table)) {
        return(NULL)
    }
    
    if (!is.null(mirna)) {
        mirna <- paste(mirna, collapse = "','")
        mirna <- paste("('", mirna, "')", sep = "")
    }
    if (!is.null(target)) {
        target <- paste(target, collapse = "','")
        target <- paste("('", target, "')", sep = "")
    }
    if (!is.null(disease.drug)) {
        disease.drug <- paste(disease.drug, collapse = "','")
        disease.drug <- paste("('", disease.drug, "')", sep = "")
    }

    if (!is.null(org)) {
        org <- tolower(org)
        if (org %in% c("hsa", "human", "homo sapiens")) {
            org <- "hsa"
        } else if (org %in% c("mmu", "mouse", "mus musculus")) {
            org <- "mmu"
        } else if (org %in% c("rno", "rat", "rattus norvegicus")) {
            org <- "rno"
        } else {
            stop(paste("Organism", org, "is not in multiMiR. Current options",
                       "are 'hsa' (human), 'mmu' (mouse) and 'rno' (rat).\n", 
                       sep = " "))
        }
    }
    
    table <- tolower(table)
    
    # prepare query for validated target table
    if (table %in% c("mirecords", "mirtarbase", "tarbase")) {
        q <- paste("SELECT m.mature_mirna_acc, m.mature_mirna_id,",
                   "t.target_symbol, t.target_entrez, t.target_ensembl,",
                   "i.experiment, i.support_type, i.pubmed_id FROM", 
                   mirna.table, "AS m INNER JOIN", table, "AS i INNER JOIN",
                   target.table, 
                   "AS t ON (m.mature_mirna_uid=i.mature_mirna_uid AND",
                   "i.target_uid=t.target_uid) WHERE", 
                   sep = " ")
        if (is.null(mirna)) {
            # only target as query
            q <- paste(q, "(t.target_symbol IN", target, 
                       "OR t.target_entrez IN", target, 
                       "OR t.target_ensembl IN", 
                       target, ")", sep = " ")
        } else if (is.null(target)) {
            # only mirna as query
            q <- paste(q, "(m.mature_mirna_acc IN", mirna, 
                       "OR m.mature_mirna_id IN", mirna, ")", 
                       sep = " ")
        } else {
            # both mirna & target as query
            q <- paste(q, "(m.mature_mirna_acc IN", mirna, 
                       "OR m.mature_mirna_id IN", mirna, 
                       ") AND (t.target_symbol IN", 
                       target, "OR t.target_entrez IN", target, 
                       "OR t.target_ensembl IN", target, ")", 
                       sep = " ")
        }
        if (!is.null(org)) {
            if (table == "tarbase") {
                q <- paste(q, " AND m.org = '", org, "'", sep = "")
            } else {
                q <- paste(q, " AND m.org = '", org, "' AND t.org = '", org,
                           "'", sep = "")
            }
        }
    }
    
    # prepare query for predicted target table
    if (table %in% c("diana_microt", "elmmo", "microcosm", "miranda", "mirdb",
                     "pictar", "pita", "targetscan")) {
        if (org == "rno" & table %in% c("diana_microt", "pictar", "pita",
                                        "targetscan")) {
            return(NULL)
        }
        q <- paste("SELECT m.mature_mirna_acc, m.mature_mirna_id,",
                   "t.target_symbol, t.target_entrez, t.target_ensembl FROM", 
                   mirna.table, "AS m INNER JOIN", table, "AS i INNER JOIN",
                   target.table, "AS t ON (m.mature_mirna_uid=i.mature_mirna_uid",
                   "AND i.target_uid=t.target_uid) WHERE", 
                   sep = " ")
        if (table == "diana_microt") {
            q <- sub(" FROM ", ", i.miTG_score AS score FROM ", q)
        } else if (table == "elmmo") {
            q <- sub(" FROM ", ", i.p AS score FROM ", q)
        } else if (table %in% c("microcosm", "mirdb", "pictar")) {
            q <- sub(" FROM ", ", i.score FROM ", q)
        } else if (table == "miranda") {
            q <- sub(" FROM ", ", i.mirsvr_score AS score FROM ", q)
        } else if (table == "pita") {
            q <- sub(" FROM ", ", i.ddG AS score FROM ", q)
        } else if (table == "targetscan") {
            q <- sub(" FROM ", ", i.context_plus_score AS score FROM ", q)
        }
        if (is.null(mirna)) {
            # only target as query
            q <- paste(q, "(t.target_symbol IN", target, 
                       "OR t.target_entrez IN", target, 
                       "OR t.target_ensembl IN", target, ")", 
                       sep = " ")
        } else if (is.null(target)) {
            # only mirna as query
            q <- paste(q, "(m.mature_mirna_acc IN", mirna, 
                       "OR m.mature_mirna_id IN", mirna, ")", 
                       sep = " ")
        } else {
            # both mirna & target as query
            q <- paste(q, "(m.mature_mirna_acc IN", mirna, 
                       "OR m.mature_mirna_id IN", mirna, 
                       ") AND (t.target_symbol IN", target, 
                       "OR t.target_entrez IN", target, 
                       "OR t.target_ensembl IN", target, ")", 
                       sep = " ")
        }
        
        # add organism to the query
        if (!is.null(org) && org %in% c("hsa","mmu","rno")) {
            q <- paste(q, " AND m.org = '", org, "' AND t.org = '", org, "'",
                       sep = "")
        }
        
        # process conserved/nonconserved sites in miranda, pita & targetscan
        # (other tables don't have conserved/nonconserved annotation)
        name <- paste(table, org, sep = ".")
        predicted.site <- match.arg(tolower(predicted.site), 
                                    c("conserved", "nonconserved", "all"))
        if (table %in% c("miranda", "pita", "targetscan")) {
            if (predicted.site == "conserved") {
                name <- paste(name, "c1", sep = ".")
                if (table == "miranda") {
                    # cut <- 0.5
                    cut <- if (org == "mmu") 0.566 else 0.57
                    q   <- paste(q, "AND i.conservation >=", cut, sep = " ")
                } else if (table == "pita") {
                    # cut <- 0.5
                    cut <- 0.9
                    q   <- paste(q, "AND i.conservation >=", cut, sep = " ")
                } else if (table == "targetscan") {
                    q   <- paste(q, "AND i.conserved_site = 'Y'", sep = " ")
                }
            } else if (predicted.site == "nonconserved") {
                name <- paste(name, "c0", sep = ".")
                if (table == "miranda") {
                    # cut <- 0.5
                    cut <- if (org == "mmu") 0.566 else 0.57
                    q   <- paste(q, "AND i.conservation <", cut, sep = " ")
                } else if (table == "pita") {
                    # cut <- 0.5
                    cut <- 0.9
                    q   <- paste(q, "AND i.conservation <", cut, sep = " ")
                } else if (table == "targetscan") {
                    q   <- paste(q, "AND i.conserved_site = 'N'", sep = " ")
                }
            }
        }

        # get dataset-specific score cutoff
        cutoffs      <- get.multimir.cutoffs()
        score.cutoff <- NA
        if (predicted.cutoff.type == "p") {
            # percent cutoff (default = 20 %)
            if (is.null(predicted.cutoff)) {
                predicted.cutoff = 20
            } else {
                if (predicted.cutoff < 1 | predicted.cutoff > 100) {
                    stop(paste("Percent predicted cutoff (predicted.cutoff)",
                               "should be between 1 and 100.\n"))
                }
            }
            score.cutoff <- 
                cutoffs[[name]][[paste(predicted.cutoff, "%", sep = "")]]
        } else if (predicted.cutoff.type == "n") {
            # number cutoff (default = 300,000)
            if (is.null(predicted.cutoff)) {
                if (cutoffs[[name]][["count"]] >= 3e+05) {
                    score.cutoff <- cutoffs[[name]][["300000"]]
                }
            } else {
                if (predicted.cutoff < 10000) {
                    message(paste("Number predicted cutoff (predicted.cutoff)",
                                  predicted.cutoff,
                                  "may be too small. A cutoff of 10000 will be",
                                  "used instead.\n", 
                                  sep = " "))
                    # predicted.cutoff = 10000
                    score.cutoff <- cutoffs[[name]][["10000"]]
                } else if (predicted.cutoff > cutoffs[[name]][["count"]]) {
                    message(paste0("Number predicted cutoff (predicted.cutoff) ",
                                   predicted.cutoff, " is larger than the total",
                                   "number of records in table ", 
                                   table, ". All records will be queried.\n"))
                } else {
                    predicted.cutoff <-
                        as.integer(as.integer(predicted.cutoff / 10000) * 10000)
                    score.cutoff <-
                        cutoffs[[name]][[as.character(predicted.cutoff)]]
                }
            }
        }
        
        # add dataset-specific cutoff to the query
        if (!is.na(score.cutoff)) {
            if (table == "diana_microt") {
                q <- paste(q, "AND i.miTG_score >=", score.cutoff, 
                           "ORDER BY i.miTG_score DESC", sep = " ")
            } else if (table == "elmmo") {
                q <- paste(q, "AND i.p >=", score.cutoff, "ORDER BY i.p DESC",
                           sep = " ")
            } else if (table %in% c("microcosm", "mirdb", "pictar")) {
                q <- paste(q, "AND i.score >=", score.cutoff, 
                           "ORDER BY i.score DESC", sep = " ")
            } else if (table == "miranda") {
                q <- paste(q, "AND i.mirsvr_score <=", score.cutoff, 
                           "ORDER BY i.mirsvr_score", sep = " ")
            } else if (table == "pita") {
                q <- paste(q, "AND i.ddG <=", score.cutoff, "ORDER BY i.ddG",
                           sep = " ")
            } else if (table == "targetscan") {
                # q <- paste(q, 'AND i.site_type == 3', sep=' ')
                q <- paste(q, "AND i.context_plus_score <=", score.cutoff,
                           "ORDER BY i.context_plus_score", 
                           sep = " ")
            }
        }
    }
    
    # prepare query for miRNA-disease/drug table
    if (table %in% c("mir2disease", "pharmaco_mir", "phenomir")) {
        cond <- 0
        if (table == "mir2disease") {
            if (is.null(mirna) & is.null(disease.drug)) {
                return(NULL)
            }
            q <- paste("SELECT m.mature_mirna_acc, m.mature_mirna_id, 'NA' AS",
                       "target_symbol, 'NA' AS target_entrez, 'NA' AS",
                       "target_ensembl, i.disease AS disease_drug,",
                       "CONCAT_WS('. ', i.year, i.title)",
                       "AS paper_pubmedID FROM", 
                       mirna.table, "AS m INNER JOIN", table, 
                       "AS i ON (m.mature_mirna_uid=i.mature_mirna_uid) WHERE", 
                       sep = " ")
            if (!is.null(disease.drug)) {
                q <- paste(q, "i.disease IN", disease.drug, sep = " ")
                cond <- cond + 1
            }
            if (!is.null(org) && org %in% c("hsa","mmu","rno") ) {
                if (cond == 0) {
                    q <- paste(q, " m.org = '", org, "'", sep = "")
                } else {
                    q <- paste(q, " AND m.org = '", org, "'", sep = "")
                }
                cond <- cond + 1
            }
        } else if (table == "pharmaco_mir") {
            q <- paste("SELECT m.mature_mirna_acc, m.mature_mirna_id,",
                       "t.target_symbol, t.target_entrez, t.target_ensembl,",
                       "i.drug AS disease_drug, i.pubmed_id AS paper_pubmedID",
                       "FROM", 
                       mirna.table, "AS m INNER JOIN", table, "AS i INNER JOIN", 
                       target.table, 
                       "AS t ON (m.mature_mirna_uid=i.mature_mirna_uid AND",
                       "i.target_uid=t.target_uid) WHERE", 
                       sep = " ")
            if (!is.null(target)) {
                q <- paste(q, "(t.target_symbol IN", target, 
                           "OR t.target_entrez IN", target, 
                           "OR t.target_ensembl IN", target, ")", 
                           sep = " ")
                cond <- cond + 1
            }
            if (!is.null(disease.drug)) {
                if (cond == 0) {
                    q <- paste(q, "i.drug IN", disease.drug, sep = " ")
                } else {
                    q <- paste(q, "AND i.drug IN", disease.drug, sep = " ")
                }
                cond <- cond + 1
            }
            if (!is.null(org) && org %in% c("hsa","mmu","rno")) {
                if (cond == 0) {
                    q <- paste(q, " m.org = '", org, "' AND t.org = '", org, "'",
                               sep = "")
                } else {
                    q <- paste(q, " AND m.org = '", org, "' AND t.org = '", org,
                               "'", sep = "")
                }
                cond <- cond + 1
            }
        } else if (table == "phenomir") {
            if (is.null(mirna) & is.null(disease.drug)) {
                return(NULL)
            }
            q <- paste("SELECT m.mature_mirna_acc, m.mature_mirna_id, 'NA' AS",
                       "target_symbol, 'NA' AS target_entrez, 'NA' AS",
                       "target_ensembl, i.disease AS disease_drug, i.pubmed_id AS",
                       "paper_pubmedID FROM", 
                       mirna.table, "AS m INNER JOIN", table, 
                       "AS i ON (m.mature_mirna_uid=i.mature_mirna_uid) WHERE", 
                       sep = " ")
            if (!is.null(disease.drug)) {
                q <- paste(q, "(i.disease IN", disease.drug, 
                           "OR i.disease_class IN", disease.drug, ")", 
                           sep = " ")
                cond <- cond + 1
            }
            if (!is.null(org) && org %in% c("hsa","mmu","rno")) {
                if (cond == 0) {
                    q <- paste(q, " m.org = '", org, "'", sep = "")
                } else {
                    q <- paste(q, " AND m.org = '", org, "'", sep = "")
                }
                cond <- cond + 1
            }
        }
        
        if (!is.null(mirna)) {
            if (cond == 0) {
                q <- paste(q, "(m.mature_mirna_acc IN", mirna, 
                           "OR m.mature_mirna_id IN", mirna, ")", 
                           sep = " ")
            } else {
                q <- paste(q, "AND (m.mature_mirna_acc IN", mirna, 
                           "OR m.mature_mirna_id IN", mirna, ")", 
                           sep = " ")
            }
            cond <- cond + 1
        }
    }
    
    # query the database
    result <- search.multimir(url = url, query = q)
    if (!is.null(result)) 
        result <- cbind(database = table, result)
    if (table %in% c("mir2disease", "pharmaco_mir", "phenomir")) 
        result <- unique(result)

    return(result)
}




