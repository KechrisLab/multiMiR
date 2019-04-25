
# when loading packages, load runs first, then attach

# Startup messages, reports default url
.onAttach <- function(libname, pkgname) {

    msg <- paste0("Welcome to multiMiR.\n\n", 
                  "multiMiR database URL has been set to the \n",
                  "default value: ", getOption("multimir.url"),"\n")
    msg <- paste0(msg,"\nDatabase Version: ", getOption("multimir.db.version"), 
                  "  Updated: ",getOption("multimir.db.updated"),"\n")

    packageStartupMessage(msg)

}

# Set default url options on load
.onLoad <- function(libname, pkgname) {

    op <- options()

    op.multimir <- list(multimir.url  = "http://multimir.org/",
                        multimir.queries = "cgi-bin/multimir_univ.pl")

    # Only set options multimir if name doesn't exist in current global options
    toset <- !(names(op.multimir) %in% names(op))
    if (any(toset)) options(op.multimir[toset])

    # Set database version options 
    vers_table <- multimir_dbInfoVersions()
    curr_vers  <- vers_table[1, ] # only choose top row (sorted in queryDBVersions)
    set_dbversion(dbversion_row = curr_vers)

    # No warning catches necessary, parse_response() takes care of failed con
    # messages
    invisible()

}

