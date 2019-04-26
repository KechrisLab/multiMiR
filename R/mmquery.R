
#' Creates all objects needed for the legacy S3 return object and the new S4
#' object.
#'
#' @return A list of data queried, summary of results, and related input
#'   parameters.
#' @return a list for packaging by \code{as.mmquery} and \code{as.mmquery_bioc}
#' @keywords internal
#' @importFrom tibble as_data_frame
#' @importFrom purrr map
#' @importFrom purrr transpose
#' @importFrom purrr reduce
#' @importFrom purrr compact
extract_mmquery <- function(outlist, org, .args, summary = FALSE, 
                            use.tibble = FALSE) {

    # outlist structured by table (list containing data, query, table name,
    # type) restructure so organized by type (predicted/validated/diseasedrug)
    outobj <- split_by(outlist, ~ .x$type) 
    outobj <- map(outobj, ~ transpose(.x)) 
    outobj <- list(data = map(outobj, ~ reduce(.x$data, rbind)),
                   queries = map(outobj, ~ reduce(.x$query, c))) 

    # Add summary table if requested
    data_wo_null <- compact(outobj$data)

    if (summary) {
        mmsum <- multimir.summary(data_wo_null)
        if (use.tibble) mmsum <- as_data_frame(mmsum)
    } else mmsum <- data.frame()

    cutoff      <- null_to_num(.args$predicted.cutoff)
    cutoff.type <- null_to_char(.args$predicted.cutoff.type)
    site        <- null_to_char(.args$predicted.site)

    rtn <- list(validated    = null_to_df(outobj$data$validated),
                predicted    = null_to_df(outobj$data$predicted),
                disease.drug = null_to_df(outobj$data$disease.drug),
                queries      = outobj$queries,
                summary      = mmsum,
                tables       = .args$table,
                org          = .args$org,
                # not including these, because they could be very large
                # mirna  = .args$mirna,
                # target = .args$target,
                # disease.drug        = .args$disease.drug,
                predicted.cutoff      = cutoff,
                predicted.cutoff.type = cutoff.type,
                predicted.site        = site)

    return(rtn)

}


#' S3 constructor and methods for object returned by \code{get_multimir()}.
#' 
#' This package's primary user-facing object. Contains the SQL statement and the
#' returned data query, as well as a summary table depending on
#' specified option. 
#' 
#' @keywords internal
#' @return An \code{mmquery} object.
as.mmquery <- function(a_list) {

    stopifnot(all(c("validated", "predicted", "disease.drug", "queries",
                    "summary", "tables", "org", "predicted.cutoff",
                    "predicted.cutoff.type", "predicted.site") %in%
                  names(a_list)))

    # Create and return s3 object
    structure(list(validated    = a_list$validated,
                   predicted    = a_list$predicted,
                   disease.drug = a_list$disease.drug,
                   queries      = a_list$queries,
                   summary      = a_list$summary),
              class  = c("mmquery"),
              tables = a_list$table,
              org    = a_list$org,
              predicted.cutoff      = a_list$predicted.cutoff,
              predicted.cutoff.type = a_list$predicted.cutoff.type,
              predicted.site        = a_list$predicted.site)

}


#' @importFrom tibble as_data_frame
#' @rdname as.mmquery
print.mmquery <- function(x) {

    cat("MultiMiR query\n")
    cat("Validated interactions:\n")
    print(as_data_frame(x$validated), n = 5)
    cat("Predicted interactions:\n")
    print(as_data_frame(x$predicted), n = 5)
    cat("Disease/Drug associations:\n")
    print(as_data_frame(x$diseasedrug), n = 5)
    cat("Summary:\n")
    print(as_data_frame(x$summary), n = 5)

}





#' S4 constructor and methods for object returned by \code{get_multimir()}.
#'
#' This package's primary user-facing object. Contains the SQL statement and the
#' returned data query, as well as a summary table depending on
#' specified option.  Note that the returned data is now contained in a single
#' dataframe. To filter to a specific type of association or interaction,
#' \code{select} on the \code{type} variable.
#'
#' @slot data A dataframe containing validated and predicted microRNA-target
#' interactions and disease/drug assocations found.
#' @slot queries A list of queries submitted to the multiMiR SQL server.
#' @slot summary A summary dataframe of the returned microRNA dataframes
#' @slot tables A character vector of the microRNA relationship types returned
#' (validated, predicted, disease.drug, or all).
#' @slot org The selected organism (hsa/human, mmu/mouse, rno/rat).
#' @slot predicted.cutoff An integer giving a prediction score cutoff.
#' @slot predicted.cutoff.type A character indicating the type of prediction
#' score cutoff (p = percentage, n = number, character() = none)
#' @slot predicted.site A character string indicating the type of predicted
#' target sites to searched.
#' @param x,object An mmquery_bioc object.
#' @param keys A result of the keys() function. For the mmquery_bioc class this
#'   is a character vector of microRNA's in the returned mmquery_bioc object.
#' @param keytype allows the user to discover which keytypes can be passes in to
#' select or keys  and the keytype argument
#' @param columns lists the columns that can be returned for the
#' \code{mmquery_bioc} object.
#' @param .list a list of returned dataframes, summary
#' @param ... additional arguments
#' @return an s4 object of class mmquery_bioc.  Contains queried data, a summary
#' dataset, and associated input parameters.
#' 
#' @importFrom AnnotationDbi select
#' @importFrom AnnotationDbi columns
#' @importFrom AnnotationDbi keys
#' @importFrom AnnotationDbi keytypes
#' @importFrom AnnotationDbi show
#' @importFrom BiocGenerics Filter
#' @importFrom methods new
#' @importFrom methods slot
#' @importFrom purrr reduce
#' @importFrom dplyr full_join
#' @export
setClass("mmquery_bioc", 
         representation(data = "data.frame",
                        queries      = "list",
                        summary      = "data.frame",
                        tables = "character",
                        org    = "character",
                        predicted.cutoff      = "numeric",
                        predicted.cutoff.type = "character",
                        predicted.site        = "character"))

#' @rdname mmquery_bioc-class
as.mmquery_bioc <- function(.list) {

    #stopifnot(all(c("validated", "predicted", "disease.drug", "queries",
    #                "summary", "table", "org", "predicted.cutoff",
    #                "predicted.cutoff.type", "predicted.site") %in%
    #              names(a_list)))
    tables <- sapply(table_types(), function(y) {
                       rtn <- .list[[y]]
                       if (nrow(rtn) > 0) rtn$type <- y
                       rtn
           })
    tables <- Filter(length, tables)
    if (length(tables) != 0) {
        tables <- reduce(tables, full_join)
    } else {
        tables <- data.frame()
    }

    # Create and return s3 object
    new("mmquery_bioc",
        data    = tables,
        queries      = .list$queries,
        summary      = .list$summary,
        tables       = .list$table,
        org          = .list$org,
        predicted.cutoff      = .list$predicted.cutoff,
        predicted.cutoff.type = .list$predicted.cutoff.type,
        predicted.site        = .list$predicted.site
        )

}

#' @rdname mmquery_bioc-class
#' @export
setMethod("columns", "mmquery_bioc", function(x) mm_cols(x))

#' @rdname mmquery_bioc-class
#' @export
setMethod("keys", "mmquery_bioc",
          function(x, keytype, ...) {

              if(missing(keytype)){
                  keytype <- mm_centralPkgSymbol()
              }
              unique(unname(x@data[[keytype]]))

          })

#' @rdname mmquery_bioc-class
#' @export
setMethod("keytypes", "mmquery_bioc", function(x) mm_cols(x))

#' @rdname mmquery_bioc-class
#' @export
setMethod("select", "mmquery_bioc",
          function(x, keys, columns, keytype, ...) {

              if (missing(keytype)) keytype <- mm_centralPkgSymbol()
              rtn <- slot(x, "data")
              rtn <- rtn[, colnames(rtn) %in% unique(c(columns, keytype))]
              if (nrow(rtn) > 0) rtn <- rtn[rtn[[keytype]] %in% keys, ]
              rtn

          })

#' @rdname mmquery_bioc-class
#' @export
setMethod("show", "mmquery_bioc",
    function(object) {

        cat("MultiMiR query\n")
        cat("Validated interactions:\n")
        print(as_data_frame(select(object, keytype = "type",
                                   keys = "validated",
                                   columns = columns(object), n = 5)))
        cat("Predicted interactions:\n")
        print(as_data_frame(select(object, keytype = "type",
                                   keys = "predicted",
                                   columns = columns(object), n = 5)))
        cat("Disease/Drug associations:\n")
        print(as_data_frame(select(object, keytype = "type",
                                   keys = "disease.drug",
                                   columns = columns(object), n = 5)))
        cat("Summary:\n")
        print(as_data_frame(object@summary), n = 5)

    })


#' @keywords internal
mm_cols <- function(x) {
    rtn <- colnames(slot(x, "data"))
    rtn <- unique(unname(rtn))
    rtn
}

#' @keywords internal
mm_centralPkgSymbol <- function() "mature_mirna_id"


# # Optional functions
# setMethod("metadata", "AnnotationDb",
#     function(x) dbReadTable(dbconn(x), "metadata")
# )
