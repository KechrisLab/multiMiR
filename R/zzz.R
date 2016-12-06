


# Startup messages, basic for now
.onAttach <- function(libname, pkgname) {
    packageStartupMessage("Welcome to multiMiR.")
}


.onLoad <- function(libname, pkgname) {
    op <- options()

    op.devtools <- list(
        multimir.url         = "http://multimir.ucdenver.edu/cgi-bin/multimir.pl",
        multimir.schema.url  = "http://multimir.ucdenver.edu/multiMiR_DB_schema.sql",
        multimir.cutoffs.url = "http://multimir.ucdenver.edu/multimir_cutoffs.rda"
        )

    # Only set options multimir if name doesn't exist in current global options
    toset <- !(names(op.devtools) %in% names(op))
    if (any(toset)) options(op.devtools[toset])

    invisible()
}
