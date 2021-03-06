\name{otosearch3}
\alias{otosearch3}
\title{Otosearch - internal function}
\usage{
otosearch3(specimen, project, show = 5)
}
\arguments{
  \item{specimen}{matrix of configuration of query}

  \item{project}{a project object, see
  \code{\link{saveproj}}}

  \item{show}{how many search results to be shown (in order
  of ranking)}
}
\value{
search result
}
\description{
search one specimen against whole database, with exhaustive
search using all possible semi-landmarks configuration
arrangements
}
\seealso{
Functions that wraps this function: \code{\link{otosearch}}

Which this function wraps: \code{\link{otosearch2}},
\code{\link{reland}}
}
\keyword{internal}

