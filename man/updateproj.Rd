\name{updateproj}
\alias{updateproj}
\title{Update the project}
\usage{
updateproj(name, landmark = NULL, des = NULL, gpa = NULL, nef = NULL,
  pc = NULL, har = NULL, des.threshold = NULL, nef.threshold = NULL,
  gpa.threshold = NULL)
}
\arguments{
  \item{name}{character. file name of the project to be
  updated.}

  \item{landmark}{p x k x n array of semi-landmarks, e.g.
  \code{A} value of \code{\link{routine2}}.}

  \item{des}{a \code{\link{routine1}} object}

  \item{gpa}{a \code{\link{rGPA}} object}

  \item{nef}{a \code{\link{rNEF}} object}

  \item{pc}{numeric. number of optimized PC range
  determined by \code{\link{pccv}}}

  \item{har}{numeric. number of optimized harmonics range
  determined by \code{\link{harcv}}}

  \item{des.threshold}{numeric. optional. threshold value
  determined by \code{\link{threcv}}}

  \item{nef.threshold}{numeric. optional. threshold value
  determined by \code{\link{threcv}}}

  \item{gpa.threshold}{numeric. optional. threshold value
  determined by \code{\link{threcv}}}
}
\value{
the selected values are updated (overwrited) with the new
addtions.
}
\description{
update the project previously saved using
\code{\link{saveproj}}
}
\details{
the objects assigned with new values will be overwrited
with new ones, the rest will remain the same.
}
\note{
the current version only support overwrite(!), not append.
}
\seealso{
\code{\link{saveproj}}, \code{\link{loadproj}}
}

