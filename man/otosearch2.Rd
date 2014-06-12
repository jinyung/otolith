\name{otosearch2}
\alias{otosearch2}
\title{Otosearch - internal}
\usage{
otosearch2(specimen, database, species = NULL, show = 5)
}
\arguments{
  \item{specimen}{matrix of configuration of query}

  \item{database}{array of configuration to search against}

  \item{species}{species/ grouping of \code{database} in
  the same order}

  \item{show}{how many search results to be shown (in order
  of ranking)}
}
\value{
search result
}
\description{
search one specimen against whole database
}
\seealso{
Functions that wraps this function:
\code{\link{otosearch3}}

Which this function wraps: \code{\link{rotonto}},
\code{\link{kendalldist}}
}
\keyword{internal}

