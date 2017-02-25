
# Startup messages, reports default url
.onAttach <- function(libname, pkgname) {

    msg <- paste0("Welcome to multiMiR.\n\n", 
                  "multiMiR database URL has been set to the \n",
                  "default value: ", getOption("multimir.url"),"\n")
    msg <- paste0(msg,"\nDatabase Version:", getOption("multimir.db.version"), 
                  "  Updated:",getOption("multimir.db.updated"),"\n")

    packageStartupMessage(msg)

}


# Set default url options on load
.onLoad <- function(libname, pkgname) {

    op <- options()
   
    op.multimir <- 
        list(
             multimir.db.version  = "0",
             multimir.db.updated  = "",
             multimir.db.name     = "",
             multimir.db.tables   = "http://multimir.ucdenver.edu/multiMiR_dbTables.txt",
             multimir.url         = "http://multimir.ucdenver.edu/cgi-bin/multimir_univ.pl",
             multimir.schema.url  = "http://multimir.ucdenver.edu/multiMiR_DB_schema.sql",
             multimir.cutoffs.url = "http://multimir.ucdenver.edu/",
             multimir.error.msg   = 
                 paste("There was an error setting up the database version.  This",
                       "indicates either a server error or that the server is currently",
                       "unavailable.  Please try again later.", sep = "\n")     
             )
    
    # Only set options multimir if name doesn't exist in current global options
    toset <- !(names(op.multimir) %in% names(op))
    if (any(toset)) options(op.multimir[toset])

    tryCatch({

        mmurl  <- "http://multimir.ucdenver.edu/cgi-bin/multimir_univ.pl"
        query  <- paste("SELECT * FROM multimir_versions.version",
                        "WHERE public=1 ORDER BY version DESC")
        result <- postForm(mmurl, query = query, .cgifields = c("query"))
        result <- readHTMLTable(result)

        current <- result[[2]][1, ]
        options(multimir.db.version  = as.character(current[[1]]))
        options(multimir.db.updated  = as.character(current[[2]]))
        options(multimir.db.name     = as.character(current[[4]]))
        options(multimir.db.tables   = paste0("http://multimir.ucdenver.edu/",
                                              as.character(current[[7]])))
        options(multimir.url         = mmurl)
        options(multimir.schema.url  = paste0("http://multimir.ucdenver.edu/",
                                              as.character(current[[5]])))
        options(multimir.cutoffs.url = paste0("http://multimir.ucdenver.edu/",
                                              as.character(current[[3]])))
        options(multimir.error.msg   = "")
    }, warning = function(war) {
        warning("There was an error setting up the database ",
                "version. This indicates either a server error or ",
                "that the server is currently unavailable. Please ",
                "try again later.")
    }, error = function(e) {
        warning("There was an error setting up the database ",
                "version. This indicates either a server error or ",
                "that the server is currently unavailable. Please ",
                "try again later.")
    }, finally = {}
    )

    invisible()

}

