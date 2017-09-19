
#' Constructor and methods for object returned by \code{get_multimir()}.
#' 
#' This package's primary user-facing object. Contains the SQL statement and the
#' returned data query, as well as a summary table depending on
#' specified option. 
#' 
#' @keywords internal
#' @return An \code{mmquery} object.
as_mmquery <- function(outlist, org, .args, summary = FALSE, 
                       use.tibble = FALSE) {

    # outlist structured by table (list containing data, query, table name,
    # type) restructure so organized by type (predicted/validated/diseasedrug)
    outobj <- split_by(outlist, ~ .x$type) 
    outobj <- purrr::map(outobj, ~ purrr::transpose(.x)) 
    outobj <- list(data = purrr::map(outobj, ~ purrr::reduce(.x$data, rbind)),
                   queries = purrr::map(outobj, ~ purrr::reduce(.x$query, c))) 

    # Add summary table if requested
    data_wo_null <- purrr::compact(outobj$data)
    if (summary) {
        mmsum <- multimir.summary(data_wo_null)
        if (use.tibble) mmsum <- tibble::as_data_frame(mmsum)
    } else mmsum <- NULL

    # Create and return s3 object
    structure(list(validated    = outobj$data$validated,
                   predicted    = outobj$data$predicted,
                   disease.drug = outobj$data$disease.drug,
                   queries      = outobj$queries,
                   summary      = mmsum),
              class  = c("mmquery"),
              tables = .args$table,
              org    = .args$org,
              # not including these, because they could be very large
              # mirna  = .args$mirna,
              # target = .args$target,
              # disease.drug          = .args$disease.drug,
              predicted.cutoff      = .args$predicted.cutoff,
              predicted.cutoff.type = .args$predicted.cutoff.type,
              predicted.site        = .args$predicted.site
              )
}

print <- function(x, ...) {
    UseMethod("print", x)
}

print.mmquery <- function(x) {

    cat("MultiMiR query\n")
    cat("Validated interactions:\n")
    tibble::as_data_frame(x$validated) 
    cat("Predicted interactions:\n")
    tibble::as_data_frame(x$predicted) 
    cat("Disease/Drug associations:\n")
    tibble::as_data_frame(x$diseasedrug) 
    cat("Summary:\n")
    tibble::as_data_frame(x$summary) 

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
# 
# sumstats <- function(x) {
#     if (is.null(x)) {
#         present <- FALSE
#         nrows   <- 0
#     } else {
#         present <- TRUE
#         nrows   <- nrow(x)
#     }
# }
