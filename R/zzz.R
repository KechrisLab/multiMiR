
# Startup messages, reports default url
.onAttach <- function(libname, pkgname) {
  msg=paste0("Welcome to multiMiR.\n\n", 
             "multiMiR database URL has been set to the ",
             "default value: ", getOption("multimir.url"),"\n")
  #if(nchar(getOption("multimir.error.msg"))>0){
  #  msg=paste0(msg,"\nError: ",getOption("multimir.error.msg"),"\n")
  #}else{
    msg=paste0(msg,"\nDatabase Version:", getOption("multimir.db.version"), 
               "  Updated:",getOption("multimir.db.updated"),"\n")
  #}
  packageStartupMessage(msg)

}


# Set default url options on load
.onLoad <- function(libname, pkgname) {

    op <- options()
   
    
    op.devtools = list(
      multimir.db.version  = "0",
      multimir.db.updated  = "",
      multimir.db.name     = "",
      multimir.db.tables  = "http://multimir.ucdenver.edu/multiMiR_dbTables.txt",
      multimir.url = "http://multimir.ucdenver.edu/cgi-bin/multimir_univ.pl",
      multimir.schema.url  = "http://multimir.ucdenver.edu/multiMiR_DB_schema.sql",
      multimir.cutoffs.url = "http://multimir.ucdenver.edu/",
      multimir.error.msg   = "The multiMiR Server did not return a result.  This is most likely from an incorrect version number."
    )
    
    # Only set options multimir if name doesn't exist in current global options
    toset <- !(names(op.devtools) %in% names(op))
    if (any(toset)) options(op.devtools[toset])
    
    tryCatch({
      mmurl="http://multimir.ucdenver.edu/cgi-bin/multimir_univ.pl"
      result <- postForm(mmurl, query = "Select * from multimir_versions.version where public=1 order by version DESC", .cgifields = c("query"))
      result <- readHTMLTable(result)
      current <- result[[2]][1,]
      options(multimir.db.version  = as.character(current[[1]]))
      options(multimir.db.updated  = as.character(current[[2]]))
      options(multimir.db.name     = as.character(current[[4]]))
      options(multimir.db.tables  = paste0("http://multimir.ucdenver.edu/",as.character(current[[7]])))
      options(multimir.url = mmurl)
      options(multimir.schema.url  = paste0("http://multimir.ucdenver.edu/",as.character(current[[5]])))
      options(multimir.cutoffs.url = paste0("http://multimir.ucdenver.edu/",as.character(current[[3]])))
      options(multimir.error.msg   = "")
      
    },warning = function(war){
      cat(war)
    },error = function(e){
      cat(e)
    },finally = {})

    invisible()

}
