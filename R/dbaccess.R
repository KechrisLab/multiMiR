

# To search the multiMiR database on the web server given a MySQL query
search.multimir <- function(url = getOption("multimir.url"), 
                            query
                            ) {
    dbName=getOption("multimir.db.name")
    result <- postForm(url, query = query, dbName=dbName, .cgifields = c("query","dbName"))
    result <- readHTMLTable(result)
    result <- parse.multimir(result)
    return(result)
}

# To switch DB version to search to the specified version if one matches
multimir_switchDBVersion <- function(url = getOption("multimir.url"),dbVer=getOption("multimir.db.version")) {
  tryCatch({
    query=paste0("Select * from multimir_versions.version where version='",dbVer,"' and public=1 order by version DESC")
    result <- postForm(url, query = query, .cgifields = c("query"))
    result <- readHTMLTable(result)
    if(as.numeric(as.character(result[[1]][[1]]))==0){
      message("ERROR:Version not set.\nVersion probably doesn't match an available version.  Please use the full version as displayed in multimir_dbInfoVersions().")
    }else{
      current <- result[[2]][1,]
      options(multimir.db.version  = as.character(current[[1]]))
      options(multimir.db.updated  = as.character(current[[2]]))
      options(multimir.db.name     = as.character(current[[4]]))
      options(multimir.db.tables  = paste0("http://multimir.ucdenver.edu/",as.character(current[[7]])))
      options(multimir.url = getOption("multimir.url"))
      options(multimir.schema.url  = paste0("http://multimir.ucdenver.edu/",as.character(current[[5]])))
      options(multimir.cutoffs.url = paste0("http://multimir.ucdenver.edu/",as.character(current[[3]])))
      options(multimir.error.msg   = "")
      cat( paste0("Now using database version: ",getOption("multimir.db.version") ))
  }}, warning = function(war){
      message(war)
  }, error = function(e){
      message(e)
  }, finally = {})
  
}

# To count records in the database
multimir_dbCount <- function(url = getOption("multimir.url")) {
    res <- search.multimir(url = url, query = "SELECT * FROM map_counts")
    for (i in 2:ncol(res)) {
        res[, i] <- as.numeric(as.character(res[, i]))
    }
    return(res)
}


# To display database information
multimir_dbInfo <- function(url = getOption("multimir.url")) {
    res <- search.multimir(url = url, query = "SELECT * FROM map_metadata")
    return(res)
}

# To display database information on DB versions
multimir_dbInfoVersions <- function(url = getOption("multimir.url")) {
  res <- search.multimir(url = url, query = "Select * from multimir_versions.version where public=1 order by version DESC")
  return(res)
}

# To display database schema
multimir_dbSchema <- function(schema.file = getOption("multimir.schema.url")) {
    schema <- readLines(schema.file)
    cat(schema, sep = "\n")
}


# To show tables in the multimir database
multimir_dbTables <- function(url = getOption("multimir.db.tables")) {
    res <- readLines(url)
    return(res)
}


# To list miRNAs, genes, drugs or diseases in the multimir database
list.multimir <- function(x   = c("mirna", "gene", "drug", "disease"),
                          url = getOption("multimir.url")) {
    x <- match.arg(tolower(x), c("mirna", "gene", "drug", "disease"))
    if (x == "mirna") {
        q <- "SELECT * FROM mirna"
    } else if (x == "gene") {
        q <- "SELECT * FROM target"
    } else if (x == "drug") {
        q <- "SELECT DISTINCT(drug) FROM pharmaco_mir"
    } else if (x == "disease") {
        # q <- 'SELECT DISTINCT(disease)FROM mir2disease UNION SELECT
        #       DISTINCT(disease) FROM phenomir' # did not work for IP's
        #       outside of campus - split the query into two
        q1 <- "SELECT DISTINCT(disease) FROM mir2disease"
        result1 <- search.multimir(url = url, query = q1)
        q2 <- "SELECT DISTINCT(disease) FROM phenomir"
        result2 <- search.multimir(url = url, query = q2)
        result  <- sort(union(toupper(result1[, 1]), toupper(result2[, 1])))
        result  <- data.frame(result)
        colnames(result) <- "disease"
        return(result)
    }
    result <- search.multimir(url = url, query = q)
    return(result)
}


# To load pre-calculated score cutoffs
get.multimir.cutoffs <- function(cutoff.file = getOption("multimir.cutoffs.url")) {
    multimir_cutoffs <- NULL
    url.file <- url(cutoff.file)
    on.exit(close(url.file))
    load(url.file)
    return(multimir_cutoffs)
}


# To parse the result from the multimir web server.  Two tables should return. The first table
# (result[[1]]) is the summary. And the second table (result[[2]]) has the result in details.
parse.multimir <- function(HTML.result) {
    result <- NULL
    l      <- length(HTML.result)
    if (l == 2) {
        result <- HTML.result[[2]]
    } else if (l == 1) {
        # cat('No records returned for your query.\n')
    } else if (l == 0) {
        cat(paste("Request to multiMiR web server failed. There could be",
                  "incorrect syntax in your query, or you are not connected",
                  "to the internet.  Alternatively the multiMiR web server",
                  "at http://multimir.ucdenver.edu is temporarily down.\n"))
    }
    return(result)
}


# The main function to search miRNA-target and miRNA-disease interactions
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


# To get miRNA-target interactions in a given table
get.multimir.by.table <- function(url = NULL,
                                  org = "hsa",
                                  mirna = NULL,
                                  target = NULL,
                                  table = NULL,
                                  disease.drug = NULL,
                                  predicted.cutoff = NULL,
                                  predicted.cutoff.type = "p",
                                  predicted.site = "conserved",
                                  mirna.table = "mirna",
                                  target.table = "target") {
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


# To summarize the result from functions get.multimir*
multimir.summary <- function(result, 
                             pair.index = 2:6, 
                             order.by = "all.sum") {
    len <- length(pair.index)
    r   <- NULL
    for (n in names(result)) {
        r <- rbind(r, cbind(result[[n]][, pair.index],
                            matrix(result[[n]]$database, ncol = 1)))
    }

    if (is.null(r)) return(NULL)
    
    info <- table(apply(r[, 1:len], 1, function(x) {
                            paste(x, collapse = "|")
                            }), r[, len + 1])
    info.ncol <- ncol(info)
    if (info.ncol > 1) {
        all.sum <- apply(info, 1, function(x) {
                             sum(x > 0)
                             })
        cols <- colnames(info)
        p.m <- match(cols, c("diana_microt", "elmmo", "microcosm", "miranda",
                             "mirdb", "pictar", "pita", "targetscan"))
        if (sum(!is.na(p.m)) > 1) {
            p.sum <- apply(matrix(info[, !is.na(p.m)], ncol = sum(!is.na(p.m))),
                           1, function(x) {
                               sum(x > 0)
                           })
            info <- cbind(info, predicted.sum = p.sum)
        } else if (sum(!is.na(p.m)) == 1) {
            p.sum <- as.integer(info[, !is.na(p.m)] > 0)
            info  <- cbind(info, predicted.sum = p.sum)
        }
        v.m <- match(cols, c("mirecords", "mirtarbase", "tarbase"))
        if (sum(!is.na(v.m)) > 1) {
            v.sum <- apply(matrix(info[, !is.na(v.m)], ncol = sum(!is.na(v.m))),
                           1, function(x) {
                               sum(x > 0)
                           })
            info <- cbind(info, validated.sum = v.sum)
        } else if (sum(!is.na(v.m)) == 1) {
            v.sum <- as.integer(info[, !is.na(v.m)] > 0)
            info  <- cbind(info, validated.sum = v.sum)
        }
        d.m <- match(cols, c("mir2disease", "pharmaco_mir", "phenomir"))
        if (sum(!is.na(d.m)) > 1) {
            d.sum <- apply(matrix(info[, !is.na(d.m)], ncol = sum(!is.na(d.m))),
                           1, function(x) {
                               sum(x > 0)
                           })
            info <- cbind(info, disease.sum = d.sum)
        } else if (sum(!is.na(d.m)) == 1) {
            d.sum <- as.integer(info[, !is.na(d.m)] > 0)
            info <- cbind(info, disease.sum = d.sum)
        }
        info <- cbind(info, all.sum = all.sum)
    }
    
    s <- NULL
    for (i in 1:nrow(info)) {
        row.name = rownames(info)[i]
        row.name = sub("\\|$", "\\|\\|", row.name)
        pair <- strsplit(row.name, "\\|")[[1]]
        pair <- c(pair, info[i, ])
        s <- rbind(s, pair)
    }
    colnames(s) <- c(colnames(result[[1]])[pair.index], colnames(info))
    s <- data.frame(s, row.names = NULL)
    
    m <- match(order.by, colnames(s))
    if (is.na(m)) {
        s <- s[order(as.numeric(as.character(s[, ncol(s)])), decreasing = TRUE), ]
    } else {
        s <- s[order(as.numeric(as.character(s[, m])), decreasing = TRUE), ]
    }
    s <- data.frame(s, row.names = NULL)
    
    n.s <- ncol(s)
    n.i <- ncol(info)
    for (n in (n.s - n.i + 1):n.s) {
        s[, n] <- as.numeric(as.character(s[, n]))
    }
    
    return(s)
}


# To change the encoding of URL (to account for the OS difference).  This is modified from the URLencode
# function in the utils package.
myurlencode <- function(url) {
    OK <- "[^-ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789$_.+*(),:/?=]"
    x <- strsplit(url, "")[[1L]]
    z <- grep(OK, x)
    if (length(z)) {
        y <- sapply(x[z], function(x) paste0("%", as.character(charToRaw(x)),
                                             collapse = ""))
        x[z] <- y
    }
    paste(x, collapse = "")
}


# To add external database link for each of the multiMiR result entry
add.multimir.links <- function(x, org) {
    warning(paste("Some of the links to external databases may be broken due",
                  "to outdated identifiers in these databases. Please refer",
                  "to Supplementary Table 2 in the multiMiR paper for details",
                  "of the issue.\n"))
    links <- rep(NA, nrow(x))
    db    <- as.character(unique(x$database))
    for (d in db) {
        m   <- which(x$database == d)
        mir <- as.character(x$mature_mirna_id[m])
        if (d == "mirecords") {
            # NOTE: need to resolve miRNA IDs with '*' in mirecords
            mir <- sub("-5p", "", mir)
            mir <- sub("-3p", "", mir)
            symbol <- as.character(x$target_symbol[m])
            if (org == "hsa") {
                s <- "species=Homo+sapiens"
            } else if (org == "mmu") {
                s <- "species=Mus+musculus"
            } else if (org == "rno") {
                s <- "species=Rattus+norvegicus"
            }
            links[m] <- paste("http://mirecords.biolead.org/interactions.php?",
                              s, "&mirna_acc=", mir,
                              "&targetgene_type=symbol&targetgene_info=", 
                              symbol, "&v=yes&search_int=Search", sep = "")
        } else if (d == "mirtarbase") {
            symbol   <- as.character(x$target_symbol[m])
            links[m] <- 
                paste("http://mirtarbase.mbc.nctu.edu.tw/php/search.php?org=",
                      org, "&mirnas=", mir, "&targets=", symbol, "&opt=adv", 
                      sep = "")
        } else if (d == "tarbase") {
            symbol   <- as.character(x$target_symbol[m])
            links[m] <- 
                paste("http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=tarbase/index&mirnas=",
                      mir, "&genes=", symbol, sep = "")
        } else if (d == "mir2disease") {
            # NOTE: Can only search by miRNA, gene or disease alone - here use gene
            symbol   <- as.character(x$target_symbol[m])
            links[m] <- 
                paste("http://watson.compbio.iupui.edu:8080/miR2Disease/searchTarget.jsp?SearchUnit=target&SearchText=",
                      symbol, "&checkbox2=Causal&checkbox2=Unspecified", 
                      sep = "")
        } else if (d == "pharmaco_mir") {
            # NOTE: Links don't work

        } else if (d == "phenomir") {
            # NOTE: search by gene
            symbol   <- as.character(x$target_symbol[m])
            links[m] <- 
                paste("http://mips.helmholtz-muenchen.de/phenomir/main/list/searchform2?query=",
                      symbol, "&selectedview=mirs&searchtype=fuzzy", 
                      sep = "")
        } else if (d == "diana_microt") {
            ensembl <- as.character(x$target_ensembl[m])
            links[m] <- 
                paste("http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=microT_CDS/results&genes=",
                      ensembl, "&mirnas=", mir, "&threshold=0", 
                      sep = "")
        } else if (d == "elmmo") {
            # NOTE: Need RefSeq accession for the gene - use miRNA only
            mir <- sub("-5p", "", mir)
            mir <- sub("-3p", "", mir)
            if (org == "hsa") {
                s <- "organism=hg"
            } else if (org == "mmu") {
                s <- "organism=mm"
            } else if (org == "rno") {
                s <- "organism=rn"
            }
            links[m] <- paste("http://www.mirz.unibas.ch/ElMMo3/?", s,
                              "&cellType=all&miRNAs[]=", mir,
                              "&predict=Predict+miRNAs+targets+!", 
                              sep = "")
        } else if (d == "microcosm") {
            mir <- sub("-5p", "", mir)
            mir <- sub("-3p", "", mir)
            if (org == "hsa") {
                s <- "genome_id=2964"
            } else if (org == "mmu") {
                s <- "genome_id=3876"
            } else if (org == "rno") {
                s <- "genome_id=5171"
            }
            symbol   <- as.character(x$target_symbol[m])
            links[m] <-
                paste("http://www.ebi.ac.uk/enright-srv/microcosm/cgi-bin/targets/v5/hit_list.pl?",
                      s, "&mirna_id=", mir, "&external_name=", symbol, 
                      sep = "")
        } else if (d == "miranda") {
            # NOTE: Could only search by gene or miRNA - use gene here
            if (org == "hsa") {
                s <- "organism=9606"
            } else if (org == "mmu") {
                s <- "organism=10090"
            } else if (org == "rno") {
                s <- "organism=10116"
            }
            symbol   <- as.character(x$target_symbol[m])
            links[m] <-
                paste("http://www.microrna.org/microrna/searchGenes.do?gene=",
                      symbol, "&", s, sep = "")
        } else if (d == "mirdb") {
            # NOTE: Could only search by gene or miRNA - use gene here
            if (org == "hsa") {
                s <- "species=Human"
            } else if (org == "mmu") {
                s <- "species=Mouse"
            } else if (org == "rno") {
                s <- "species=Rat"
            }
            symbol <- as.character(x$target_symbol[m])
            links[m] <- paste("http://mirdb.org/cgi-bin/search.cgi?", s,
                              "&searchType=gene&geneChoice=symbol&searchBox=", 
                              symbol, sep = "")
        } else if (d == "pictar") {
            # NOTE: Links don't work

        } else if (d == "pita") {
            mir <- sub("-5p", "", mir)
            mir <- sub("-3p", "", mir)
            if (org == "hsa") {
                s <- "Organism=Human"
            } else if (org == "mmu") {
                s <- "Organism=Mouse"
            }
            symbol <- as.character(x$target_symbol[m])
            links[m] <-
                paste("http://genie.weizmann.ac.il/cgi-bin/search_mir07_prediction.pl?",
                      s, "&microRNAs=", mir, "&Genes=", symbol,
                      "&MinimumSeed=7&AllowSingleGU=1&AllowSingleMismatch=1&MinConservation=0&FlankOption=0_0",
                      sep = "")
        } else if (d == "targetscan") {
            mir <- sub("-5p", "", mir)
            mir <- sub("-3p", "", mir)
            symbol <- as.character(x$target_symbol[m])
            if (org == "hsa") {
                links[m] <-
                    paste("http://www.targetscan.org/cgi-bin/targetscan/vert_61/targetscan.cgi?species=Human&gid=",
                          symbol, "&mirg=", mir, sep = "")
            } else if (org == "mmu") {
                links[m] <-
                    paste("http://www.targetscan.org/cgi-bin/targetscan/mmu_61/targetscan.cgi?species=Mouse&gid=",
                          symbol, "&mirg=", mir, sep = "")
            }
        }
    }

    x = data.frame(x, DB.link = links)
    return(x)
}

