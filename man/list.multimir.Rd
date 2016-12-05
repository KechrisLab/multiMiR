\name{list.multimir}

\alias{list.multimir}

\title{
  List microRNAs, Genes, Drugs Or Diseases in the multiMiR Package
}

\description{
  \code{list.multimir} lists all the unique microRNAs, target genes,
  drugs, or diseases in the web server of the multiMiR package.
}

\usage{
  list.multimir(x = c("mirna","gene","drug","disease"),
                url = multimir_URL)
}

\arguments{
  \item{x}{
    a character string indicating what to list. This must be one of
    the strings \code{"mirna"} (default), \code{"gene"}, \code{"drug"},
    or \code{"disease"}. This can be abbreviated and is case insensitive.
  }
  \item{url}{
    a character string for the URL of the multiMiR web server.
    The default is multimir_URL
    ("http://multimir.ucdenver.edu/cgi-bin/multimir.pl").
  }
}

\details{
  \code{list.multimir} returns all the unique microRNAs, target genes,
  drugs, or diseases in the web server of the multiMiR package.
}

\value{
  \code{list.multimir} returns a data frame with information of
  microRNAs (microRNA unique ID, organism, mature microRNA accession
  number, and mature microRNA ID), target genes (gene unique ID,
  organism, gene symbol, Entrez gene ID, and Ensembl gene ID), drugs
  (drug names), and diseases (disease name).
}

\author{
  Yuanbin Ru \email{ruyuanbin@gmail.com}
}

\examples{
  miRNAs <- list.multimir("mirna")
  genes <- list.multimir("gene")
  drugs <- list.multimir("drug")
  diseases <- list.multimir("disease")
}

\keyword{utilities}
\keyword{database}

