\name{otopred}
\alias{otopred}
\title{Predict new samples based on classification models}
\usage{
otopred(query, project, type = c("nef", "gpa", "des"), method = c("lda",
  "tree", "plsda"), har = NULL, pc = NULL, threshold = NULL,
  reland = TRUE, tol = 0.1, write = FALSE)
}
\arguments{
  \item{query}{dataframe of shape descriptor data or p x k
  x n array of semi-landmarks}

  \item{project}{the project data, see
  \code{\link{saveproj}}}

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

  \item{write}{logical. whether to save the result.}
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
