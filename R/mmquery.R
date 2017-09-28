
#' Creates all objects needed for the legacy S3 return object and the new S4
#' object.
#'
#' @return A list of data queried, summary of results, and related input
#'   parameters.
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

null_to_df <- function(x) if (is.null(x)) data.frame() else x
null_to_num <- function(x) if (is.null(x)) numeric() else x
null_to_char <- function(x) if (is.null(x)) character() else x


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
                    "summary", "table", "org", "predicted.cutoff",
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


# summary.mmquery <- function(x) {
# 
#     has_validated    <- !is.null(x$validated)
#     has_predicted    <- !is.null(x$predicted)
#     has_disease.drug <- !is.null(x$disease.drug)
# 
#     cat("MultiMiR Query of microRNA\n")
# #     cat(sprintf("%s validated records found\n", nrow_validated))
# #     cat(sprintf("%s predicted records found\n", nrow_predcted))
# #     cat(sprintf("%s disease/drug records found\n", nrow_disease))
# 
# }


# sumstats <- function(x) {
#     if (is.null(x)) {
#         present <- FALSE
#         nrows   <- 0
#     } else {
#         present <- TRUE
#         nrows   <- nrow(x)
#     }
# }




#' S4 constructor and methods for object returned by \code{get_multimir()}.
#' 
#' This package's primary user-facing object. Contains the SQL statement and the
#' returned data query, as well as a summary table depending on
#' specified option. 
#' 
#' @importFrom methods new
#' @keywords internal
#' @return An \code{mmquery_bioc} object.
#' 
as.mmquery_bioc <- function(a_list) {

    #stopifnot(all(c("validated", "predicted", "disease.drug", "queries",
    #                "summary", "table", "org", "predicted.cutoff",
    #                "predicted.cutoff.type", "predicted.site") %in%
    #              names(a_list)))

    # Create and return s3 object
    new("mmquery_bioc",
        validated    = a_list$validated,
        predicted    = a_list$predicted,
        disease.drug = a_list$disease.drug,
        queries      = a_list$queries,
        summary      = a_list$summary,
        tables       = a_list$table,
        org          = a_list$org,
        predicted.cutoff      = a_list$predicted.cutoff,
        predicted.cutoff.type = a_list$predicted.cutoff.type,
        predicted.site        = a_list$predicted.site
        )

}

setClass("mmquery_bioc", 
         representation(validated    = "data.frame",
                        predicted    = "data.frame",
                        disease.drug = "data.frame",
                        queries      = "list",
                        summary      = "data.frame",
                        tables = "character",
                        org    = "character",
                        predicted.cutoff      = "numeric",
                        predicted.cutoff.type = "character",
                        predicted.site        = "character"))


#' S4 methods for mmquery_bioc class - based on AnnotationDbi accessor methods
#'
#' @param x An mmquery_bioc object.
#' @param keys A result of the keys() function. For the mmquery_bioc class this
#'   is a character vector of microRNA's in the returned mmquery_bioc object.
#' @param keytype allows the user to discover which keytypes can be passes in to
#' select or keys  and the keytype argument
#' @param columns lists the columns that can be returned for the
#' \code{mmquery_bioc} object.
#' @param ... additional arguments
#'
#' @importFrom AnnotationDbi columns
#' @importFrom AnnotationDbi keys
#' @importFrom methods slot
#' @export
setMethod("columns", "mmquery_bioc",
          function(x) {
              tables <- c("validated", "predicted", "disease.drug")
              rtn <- sapply(tables, function(y) colnames(slot(x, y)))
              rtn <- unique(unname(unlist(rtn)))
              rtn
          })

#' @rdname columns
#' @export
# keys likely miRNA
setMethod("keys", "mmquery_bioc",
          function(x, keytype, ...) {
              if(missing(keytype)){
                  keytype <- "mature_mirna_id"
              }
              tables <- c("validated", "predicted", "disease.drug")
              rtn <- sapply(tables, function(y) {
                                tbl <- slot(x, y)
                                tbl[[keytype]]
                        })
              unique(unname(unlist(rtn)))

          })


# setMethod("select", "mmquery_bioc",
#           function(x, keys, columns, keytype, ...) {
#               if (missing(keytype)) keytype <- chooseCentralOrgPkgSymbol(x)
#               jointype <- .chooseJoinType(x)
#               .select(x, keys, columns, keytype, jointype=jointype, ...)
#               ## .selectWarnJT(x, keys, columns, keytype, jointype=jointype,
#               ##               kt=kt, ...)
#           })

# #' @rdname select
# 
# #' @rdname select
# setMethod("keytypes", "mmquery_bioc",
#     function(x) {
#         kts <- .cols(x, baseType="ENTREZID")
#         .filterDeprecatedKeytypes(kts)
#     })
# 
# #' @rdname select
# setMethod("select", "mmquery_bioc",
#           function(x, keys, columns, keytype, ...) {
#               if (missing(keytype)) keytype <- chooseCentralOrgPkgSymbol(x)
#               jointype <- .chooseJoinType(x)
#               .select(x, keys, columns, keytype, jointype=jointype, ...)
#               ## .selectWarnJT(x, keys, columns, keytype, jointype=jointype,
#               ##               kt=kt, ...)
#           })
# 
# setMethod("show", "mmquery_bioc",
#     function(object)
#     {
#         cat(class(object), "object:\n")
#         metadata <- metadata(object)
#         for (i in seq_len(nrow(metadata))) {
#             cat("| ", metadata[i, "name"], ": ", metadata[i, "value"],
#                 "\n", sep="")
#         }
#         message("\n","Please see: help('select') for usage information", sep="")
#     }
# )
# 
# # Optional functions
# setMethod("metadata", "AnnotationDb",
#     function(x) dbReadTable(dbconn(x), "metadata")
# )
