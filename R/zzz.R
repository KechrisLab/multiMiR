
# Startup messages, reports default url
.onAttach <- function(libname, pkgname) {
  
    packageStartupMessage(paste0("Welcome to multiMiR.\n\n", 
                                 "multiMiR database URL has been set to the ",
                                 "default value: ", getOption("multimir.url")))

}


# Set default url options on load
.onLoad <- function(libname, pkgname) {

    op <- options()

    op.devtools <- list(
        multimir.db.version  = "2.0",
        multimir.url         = "http://multimir.ucdenver.edu/cgi-bin/multimir.pl",
        multimir.schema.url  = "http://multimir.ucdenver.edu/multiMiR_DB_schema.sql",
        multimir.cutoffs.url = "http://multimir.ucdenver.edu/multimir_cutoffs.rda"
        )

    # Only set options multimir if name doesn't exist in current global options
    toset <- !(names(op.devtools) %in% names(op))
    if (any(toset)) options(op.devtools[toset])

    invisible()

}
