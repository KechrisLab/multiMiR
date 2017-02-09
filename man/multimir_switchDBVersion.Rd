\name{multimir_switchDBVersion}

\alias{multimir_switchDBVersion}
\alias{multimir_dbInfoVersions}


\title{
  Manage Database Version to use
}

\description{
  Functions for managing the database version used to complete requests on the web server.
}

\usage{
  multimir_dbInfoVersions(url = getOption("multimir.url"))
  multimir_switchDBVersion(url = getOption("multimir.url"),dbVer="2.1.0")
}

\arguments{
  \item{url}{
    A character string naming the URL of the web server hosting the multiMiR
    database.
  }
  \item{dbVer}{
    A character string containing the full version number for the database version to 
    use for for all package functions. The default will be the most recent version.
  }
}

\details{
  \code{url} is a character string containing the URL of the multiMiR
  web server. Optional as it is set when the package is loaded.

  \code{multimir_dbInfoVersions} returns other information about the multiMiR
  database versions available. This provides a list of available options if 
  switching to previous version is desired.
  
  \code{multimir_switchDBVersion} returns other information about the multiMiR
  database versions available. This provides a list of available options if 
  switching to previous version is desired.
  
}

\value{
  
  \code{multimir_dbInfoVersions}: a data frame with information about the multiMiR
  database versions.

  \code{multimir_switchDBVersion}: none (invisible \code{NULL}).

}

\author{
  Yuanbin Ru \email{ruyuanbin@gmail.com}
}

\examples{
  multimir_dbInfoVersions()
    VERSION    UPDATED    RDA                      DBNAME          SCHEMA                 PUBLIC   TABLES
1   2.1.0      2016-12-22 multimir_cutoffs_2.1.rda multimir2_1     multiMiR_DB_schema.sql      1   multiMiR_dbTables.txt
2   2.0.0      2015-05-01 multimir_cutoffs.rda     multimir        multiMiR_DB_schema.sql      1   multiMiR_dbTables.txt
  
  multimir_switchDBVersion(dbVer="2.0.0")
  Now using database version: 2.0.0
  
}

\keyword{utilities}
\keyword{database}

