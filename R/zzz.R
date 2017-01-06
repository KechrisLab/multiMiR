
# Startup messages, reports default url
.onAttach <- function(libname, pkgname) {
  msg=paste0("Welcome to multiMiR.\n\n", 
             "multiMiR database URL has been set to the ",
             "default value: ", getOption("multimir.url"),"\n")
  if(nchar(getOption("multimir.error.msg"))>0){
    msg=paste0(msg,"\nError: ",getOption("multimir.error.msg"),"\n")
  }else{
    msg=paste0(msg,"\nDatabase Version:", getOption("multimir.db.version"), 
               "  Updated:",getOption("multimir.db.updated"),"\n")
  }
  packageStartupMessage(msg)

}


# Set default url options on load
.onLoad <- function(libname, pkgname) {

    op <- options()
   
    
    op.devtools = tryCatch({
      mmurl="http://multimir.ucdenver.edu/cgi-bin/multimir_univ.pl"
      result <- postForm(mmurl, query = "Select * from multimir_versions.version order by version DESC", .cgifields = c("query"))
      result <- readHTMLTable(result)
      current <- result[[2]][1,]
      tmp<-list(
        multimir.db.version  = as.character(current[[1]]),
        multimir.db.updated  = as.character(current[[2]]),
        multimir.url = mmurl,
        multimir.schema.url  = paste0("http://multimir.ucdenver.edu/",as.character(current[[5]])),
        multimir.cutoffs.url = paste0("http://multimir.ucdenver.edu/",as.character(current[[3]])),
        multimir.db.name     = as.character(current[[4]]),
        multimir.error.msg   = ""
      )
    },warning = function(war){
      print(war)
      return(list(
        multimir.db.version  = "0",
        multimir.db.updated  = "",
        multimir.url = mmurl,
        multimir.schema.url  = "http://multimir.ucdenver.edu/multiMiR_DB_schema.sql",
        multimir.cutoffs.url = "http://multimir.ucdenver.edu/",
        multimir.db.name     = "",
        multimir.error.msg   = "The multiMiR Server did not respond with a list of versions.  The server is temporarily unavailable.  Please try again later."
      ))
    },error = function(e){
      print(e)
      return(list(
        multimir.db.version  = "0",
        multimir.db.updated  = "",
        multimir.url = mmurl,
        multimir.schema.url  = "http://multimir.ucdenver.edu/multiMiR_DB_schema.sql",
        multimir.cutoffs.url = "http://multimir.ucdenver.edu/",
        multimir.db.name     = "",
        multimir.error.msg   = "The multiMiR Server did not respond with a list of versions.  The server is temporarily unavailable.  Please try again later."
      ))
    },finally = {})
   
    
    # Only set options multimir if name doesn't exist in current global options
    toset <- !(names(op.devtools) %in% names(op))
    if (any(toset)) options(op.devtools[toset])

    invisible()

}
