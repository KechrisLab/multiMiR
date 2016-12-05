\name{multimir_dbInfo}

\alias{multimir_URL}
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
  multimir_URL
  multimir_dbCount(url = multimir_URL)
  multimir_dbInfo(url = multimir_URL)
  multimir_dbSchema(schema.file = "http://multimir.ucdenver.edu/multiMiR_DB_schema.sql")
  multimir_dbTables(url = multimir_URL)
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
  \code{multimir_URL} is a character string containing the URL of the multiMiR
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
  \code{multimir_URL}: a character string with the URL of the multiMiR
  web server.

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
  multimir_URL

  db.count <- multimir_dbCount()

  db.info <- multimir_dbInfo()

  multimir_dbSchema()

  db.tables <- multimir_dbTables()
}

\keyword{utilities}
\keyword{database}

