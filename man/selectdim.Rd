\name{selectdim}
\alias{selectdim}
\title{Select training data dimension}
\usage{
selectdim(train, har, pc)
}
\arguments{
  \item{train}{Matrix of EFA coefficients or PC scores
  result from \code{\link{rNEF}} or \code{\link{rGPA}}
  function (called as \code{object$coeff} or
  \code{object$score}) to be trimmed down}

  \item{har}{numeric. Number of harmonics (for coefficient
  matrix)}

  \item{pc}{numeric. Number of pc (for PC score matrix)}
}
\value{
coefficients/ PC scores of selected number ready to be used
for training
}
\description{
a utility function to trim down the matrix of coefficients
(NEF)/ PC scores (GPA)
}
\seealso{
input from: \code{\link{rNEF}}, \code{\link{rGPA}}

optimising \code{pc/har} number: \code{\link{pccv}},
\code{\link{harcv}}
}

