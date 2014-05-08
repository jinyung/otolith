\name{otosearch}
\alias{otosearch}
\title{Search new specimens' configuration against project database}
\usage{
otosearch(query, project, show = 5, save = FALSE, name = NULL)
}
\arguments{
  \item{query}{p x k matrix or p x k x n array. the query
  configuration(s)}

  \item{project}{a project saved using
  \code{\link{saveproj}}. query is searched against the
  whole array of configurations saved within the project.}

  \item{show}{numeric. how many search results to be shown
  (in order of ranking).}

  \item{save}{logical. whether to save the result.}

  \item{name}{character. optional. the file name, if
  \code{save=TRUE}.}
}
\value{
The search result. if \code{save=TRUE}, the result is
written into a \code{.txt} file
}
\description{
search the semi-landmarks configuration(s) of new, unknown
specimens against the configurations saved in the project.
}
\details{
the search is based on the Procrustes (Riemannian) distance
between the query and the database. the lower the distance
(\code{rdist}) between the query and the project's
configuration, the higher their ranking in the search
result. Perfect match will have 0 distance.

To shorten the time of searching, the query is first
searched against the meanshape configurations of each
species present in the project. For each species in the
project, maximum distance between individuals and the
species's meanshape was already calculated (by
\code{\link{sprdist}}, wrapped within \code{\link{rGPA}}).
First the search will determine if the query to each
species's meanshape distance is within this range, if say,
the query is within the range of 2 species' range, then
searching will be continued with the individuals of these 2
species only. Otherwise, searching will be continued with
the individuals of the top 5 closest species only.

The search has taken the semi-landmarks configuration
arrangements into consideration already, hence the user
doesn't need to know the side/ direction of the query.
}
\seealso{
Similar: \code{\link{otopred}}

Which this function wraps: \code{\link{otosearch3}}
}

