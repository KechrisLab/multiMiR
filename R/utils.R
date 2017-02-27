



#' Load Pre-calculated Prediction Score Cutoffs in the multiMiR Package
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please set prediction score cutoff in \code{get.multimir}.
#' 
#' 
#' @export get.multimir.cutoffs
get.multimir.cutoffs <- function(cutoff.file = getOption("multimir.cutoffs.url")) {
    # To load pre-calculated score cutoffs
    # NOTE: should this fn not be exported?
    multimir_cutoffs <- NULL
    url.file         <- url(cutoff.file)
    on.exit(close(url.file))
    load(url.file)
    return(multimir_cutoffs)
}



#' Encode a URL Before Submitting It to the multiMiR Web Server
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get.multimir}.
#' 
#' 
#' @export myurlencode
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



#' Encode a URL Before Submitting It to the multiMiR Web Server
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get.multimir}.
#' 
#' 
#' @export myurlencode
myurlencode <- function(url) {
    # To change the encoding of URL (to account for the OS difference).  This is modified from the URLencode
    # function in the utils package.

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




#' Summarize microRNA/target Information from the multiMiR Package
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get.multimir}.
#' 
#' 
#' @export multimir.summary
multimir.summary <- function(result, 
                             pair.index = 2:6, 
                             order.by = "all.sum") {
    # To summarize the result from functions get.multimir*
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




# 
deprecate_arg <- function(name = c("url", "schema.file", "db.tables")) {

    name <- match.arg(name)
    ops  <- switch(name,
                   url         = "multimir.url",
                   schema.file = "multimir.schema.url",
                   db.tables   = "multimir.db.tables")

    # the function using the schema option had an arg name of url, so switch for
    # an accurate message
    if (name == "db.tables") name <- "url"

    message("The ", name, " argument is deprecated. Please set using the package ",
            "option ", ops, " via options()")

}
