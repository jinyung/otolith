\name{otosearch}
\alias{otosearch}
\title{Search new specimens against project database}
\usage{
otosearch(project, query, show = 5, saveresult = FALSE, showplot = TRUE,
  savename)
}
\arguments{
  \item{project}{path to a project (\code{.rds} file) to be
  read, saved using \code{\link{saveproj}}, or a project
  object already read into R. if not given, interactive
  file selector will pop out to prompt user to select a
  \code{.rds} file (Windows only)}

  \item{query}{path(s) to otolith images/ path(s) of folder
  containing the otolith images/ \code{.tps} file
  containing the semi-landmark configurations/ p x k matrix
  or p x k x n array of semi-landmark configuration(s) to
  be searched. If none is given, interactive file selector
  will pop out to prompt user to select images to be
  searched (Windows only)}

  \item{show}{integer. how many search results to be shown
  (in order of ranking). Note that if number of specimens
  in database searched against is less than show, the
  result shown will be less than the given number}

  \item{saveresult}{logical. whether to save the result}

  \item{showplot}{logical. whether to show plot of query
  and their matches}

  \item{savename}{character. optional. the file name, if
  \code{saveresult=TRUE}}
}
\value{
The search result. if \code{saveresult=TRUE}, the result is
written into a \code{.txt} file. The results include:
\itemize{ \item \code{rdist}: distance between the query
and the match \item \code{label}: label of the match from
the database \item \code{species}: species of the match
\item \code{inside}: whether the \code{rdist} is within the
range of that species in the database \item \code{orient}:
type of configuration that match }
}
\description{
search the images/ semi-landmarks configuration(s) of new,
unknown specimens against the configurations saved in the
project.
}
\details{
The search is based on the Procrustes (Riemannian) distance
between the query and the database. the lower the distance
(\code{rdist}) between the query and the project's
configuration, the higher their ranking in the search
result. Perfect match will have 0 distance (for non-sliding
landmark method).

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

\bold{Note:} sliding semi-landmark is currently not
supported due to the difficulty to achieve short searching
time and correct configuration guessing at the same time.
Therefore the \code{gpa} object saved in the project should
be built with fixed semi-landmarks. However, this does not
affect \code{\link{otopred}}, thus if sliding is preferred,
\code{otopred} will still give good predictions.
}
\seealso{
Similar: \code{\link{otopred}}
}

