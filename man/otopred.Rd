\name{otopred}
\alias{otopred}
\title{Predict new samples based on classification models}
\usage{
otopred(project, query, type = c("nef", "gpa", "des"), method = c("lda",
  "tree", "plsda"), har, pc, threshold, reland = TRUE, tol = 0.1,
  mode = c("pred", "search+pred"), write = FALSE, ...)
}
\arguments{
  \item{project}{path to a project (\code{.rds} file) to be
  read, saved using \code{\link{saveproj}}, or a project
  object already read into R.  if not given, interactive
  file selector will pop out to prompt user to select a
  \code{.rds} file (Windows only)}

  \item{query}{path(s) to otolith images/ path(s) of folder
  containing the otolith images/ code{.tps} file containing
  the semi-landmark configurations/ p x k matrix or p x k x
  n array of semi-landmark configuration(s) to be
  predicted. If none is given, interactive file selector
  will pop out to prompt user to select images to be
  searched (Windows only)}

  \item{type}{type of data to predict.}

  \item{method}{classification method. see note.}

  \item{har}{numeric. optional. By default \code{har} range
  saved in the project is used. a different value could be
  set using this argument.}

  \item{pc}{numeric. optional. By default \code{pc} range
  saved in the project is used. a different value could be
  set using this argument.}

  \item{threshold}{numeric. optional. threshold value on
  posterior probalility to reject the prediction. see
  \code{\link{threcv}}. Currently works with \code{lda}
  only.}

  \item{reland}{logical. whether to do automatic
  re-arrangement of landmark- configuration.}

  \item{tol}{numeric. max limit of distance (see
  \code{\link{sprdist}}) to which the automatic
  re-arrangement of landmark configuration should be
  carried out.}

  \item{mode}{when \code{="search+pred"}, searching is
  included in addition to prediction using
  \code{\link{otosearch}}. automatically changed to
  \code{"search+pred"} when \code{reland=TRUE}}

  \item{write}{logical. whether to save the result.}

  \item{...}{other arguments that can be paseed to
  \code{\link{otosearch}}}
}
\value{
matrix of prediction class and posterior probablity.
}
\description{
Predict new, unknown samples using either shape decriptors
or semi-landmark configurations.
}
\details{
The functions include the \code{\link{otosearch}} algorithm
as a mean to guess the semi-landmarks arragement (the four
types of arrangements, see \code{\link{reland}}), so that
the user can predict the new, unknown samples even if the
side and the direction of the otoliths are unknown. This
could be turned off using \code{reland=FALSE} to speed up
the prediction if the side and direction of the query is
known, and is already in the same arrangements as the
dataset in project.
}
\note{
Currently the \code{method} supported is limited to
\code{"lda"} only.
}
\seealso{
Which this function wraps: \code{\link{otosearch}}

Methods of classifier: \code{\link{lda}},
\code{\link{plsda}}, \code{\link{tree}}
}

