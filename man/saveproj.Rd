\name{saveproj}
\alias{saveproj}
\title{Save a project}
\usage{
saveproj(name, landmark = NULL, des = NULL, gpa = NULL, nef = NULL,
  pc = NULL, har = NULL, des.threshold = NULL, nef.threshold = NULL,
  gpa.threshold = NULL)
}
\arguments{
  \item{name}{character. file name of the project, should
  end with \code{".rds"}}

  \item{landmark}{p x k x n array of semi-landmarks, e.g.
  \code{A} value of \code{\link{routine2}}.}

  \item{des}{a \code{\link{routine1}} object}

  \item{nef}{a \code{\link{rNEF}} object}

  \item{gpa}{a \code{\link{rGPA}} object}

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
a \code{.rds} file saved
}
\description{
this function save the processed semi-landmarks, shape
descriptor data, GPA/ EFA transformed data, optimized
parameters into a \code{.rds} file.
}
\seealso{
Similar: \code{\link{readRDS}}, \code{\link{updateproj}}

Which this function wraps: \code{\link{saveRDS}}
}

