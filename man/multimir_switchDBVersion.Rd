\name{multimir_switchDBVersion}

\alias{multimir_switchDBVersion}


\title{
  Manage Database Version to use
}

\description{
  Functions for managing the database version used to complete requests on the web server.
}

\usage{
  multimir_switchDBVersion(url = getOption("multimir.url"),dbVer=getOption("multimir.db.version"))
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
  multimir_switchDBVersion(dbVer="2.0.0")
}

\keyword{utilities}
\keyword{database}

