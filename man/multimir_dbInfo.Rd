\name{multimir_dbInfo}

\alias{multimir_dbCount}
\alias{multimir_dbInfo}
\alias{multimir_dbSchema}
\alias{multimir_dbTables}

\title{
  Collect Information About the Web Server And Database of the multiMiR Package
}

\description{
  Functions for collecting and displaying information about the web server
  and database of the multiMiR package.
}

\usage{
  multimir_dbCount(url = getOption("multimir.url"))
  multimir_dbInfo(url = getOption("multimir.url"))
  multimir_dbSchema(schema.file = getOption("multimir.schema.url"))
  multimir_dbTables(url = getOption("multimir.url"))
}

\arguments{
  \item{url}{
    A character string naming the URL of the web server hosting the multiMiR
    database.
  }
  \item{schema.file}{
    A character string naming the file containing the multiMiR database schema.
  }
}

\details{
  \code{url} is a character string containing the URL of the multiMiR
  web server. 

  \code{multimir_dbCount} returns counts of records in the tables in the
  multiMiR database. Each table contains data from an external miRNA/target
  database.

  \code{multimir_dbInfo} returns other information about the multiMiR
  database. This includes information of external miRNA/target databases
  in multiMiR.

  \code{multimir_dbSchema} prints the schema definition of the multiMiR
  database.

  \code{multimir_dbTables} returns the list of tables in the multiMiR database.
}

\value{
  \code{multimir_dbCount}: a data frame with the count of records in each of
  the tables in the multiMiR database.

  \code{multimir_dbInfo}: a data frame with information about the multiMiR
  database.

  \code{multimir_dbSchema}: none (invisible \code{NULL}).

  \code{multimir_dbTables}: a data frame with table names in the multiMiR
  database.
}

\author{
  Yuanbin Ru \email{ruyuanbin@gmail.com}
}

\examples{
  getOption("multimir.url")
  options(multimir.url = "https://multimir.ucdenver.edu/")

  db.count <- multimir_dbCount()

  db.info <- multimir_dbInfo()

  multimir_dbSchema()

  db.tables <- multimir_dbTables()
}

\keyword{utilities}
\keyword{database}

