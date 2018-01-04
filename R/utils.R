
#' Load Pre-calculated Prediction Score Cutoffs in the multiMiR Package
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please set prediction score cutoff in \code{get_multimir}.
#' 
#' @param cutoff.file Deprecated. Set path to cutoffs file with the global
#' option \code{multimir.cutoffs}.
#' @return Cutoff values object from remote database.
#' @keywords internal
get.multimir.cutoffs <- function(name = NULL, cutoff.file = NULL) {
    # To load pre-calculated score cutoffs
    # NOTE: should this fn be exported? (NO)

    if (!is.null(cutoff.file)) deprecate_arg("cutoff.file")

    multimir_cutoffs <- NULL
    url.file         <- url(full_url("multimir.cutoffs"))
    on.exit(close(url.file))
    load(url.file)

    if (is.null(name)) {
        return(multimir_cutoffs)
    } else {
        return(multimir_cutoffs[[name]])
    }

}


#' Internal function for sending deprecation messages
#'
#' @param name Name of a deprecated function argument.
#' @return A message indicating deprecated arg and new version.
#' @keywords internal
deprecate_arg <- function(name = c("url", "schema.file", "db.tables",
                                   "cutoff.file")) {

    name <- match.arg(name)
    ops  <- switch(name,
                   url         = "multimir.url",
                   schema.file = "multimir.schema",
                   db.tables   = "multimir.db.tables",
                   cutoff.file = "multimir.cutoffs")

    # the function using the schema option had an arg name of url, so switch for
    # an accurate message
    if (name == "db.tables") name <- "url"

    message("The ", name, " argument is deprecated. Please set using the ",
            "package option ", ops, " via options()")

}


#' Internal function for adding single quotes around a string
#'
#' @param x a string to be wrapped in single quotes.
#' @return The input wrapped in single quotes.
#' @keywords internal
quote_wrap <- function(x) paste0("'", x, "'")



#' Prep certain names for use in SQL query by adding parens
#'
#' @return The input value wrapped in quotes and then parentheses.
#' @keywords internal
parens_quote <- function(x) {
    if (!is.null(x)) parens_wrap(quote_wrap(x))
}

#' Collapse a vector to a single comma-separated string and wrap in parentheses 
#'
#' @return The input vector converted to a comma-separated string wrapped in
#' parentheses.
#' @keywords internal
parens_wrap <- function(x) {
    paste0("(", paste(x, collapse = ", "), ")")
}

#' Pad single space on each side of an input 
#'
#' @return Input value wrapped in single spaces.
#' @keywords internal
pad <- function(x) paste0(" ", x, " ") 


#' Split, order and sort lists by their components.
#'
#' Copied from purrr:v0.2.2
#'
#' @param .x A list or atomic vector.
#' @param .f A function, formula, or atomic vector.
#' @param ... Additional arguments passed on to \code{.f}.
#' @return A list split by \code{.f}
#' @importFrom purrr map
#' @importFrom purrr simplify_all
#' @importFrom purrr transpose
#' @keywords internal
split_by <- function(.x, .f, ...) {
    vals <- purrr::map(.x, .f, ...)
    split(.x, purrr::simplify_all(purrr::transpose(vals)))
}

#' Replace nulls with an empty object of each type
#'
#' @return an empty \code{data.frame}, \code{numeric}, or \code{character}
#' vector.
#' @param x input object
#' @keywords internal
null_to_df <- function(x) if (is.null(x)) data.frame() else x

#' @rdname null_to_df
null_to_num <- function(x) if (is.null(x)) numeric() else x

#' @rdname null_to_df
null_to_char <- function(x) if (is.null(x)) character() else x

