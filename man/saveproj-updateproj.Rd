\name{saveproj / updateproj}
\alias{saveproj}
\alias{updateproj}
\title{Save / update project}
\usage{
saveproj(name, landmark = NULL, des = NULL, gpa = NULL, nef = NULL,
  pc = NULL, har = NULL, class = NULL)
updateproj(name, landmark = NULL, des = NULL, gpa = NULL, nef = NULL,
  pc = NULL, har = NULL, class = NULL)
}
\arguments{
  \item{name}{character or path. file name/ path of the project to be saved or to 
  be updated, should end with \code{".rds"}}

  \item{landmark}{p x k x n array of semi-landmarks, e.g.
  \code{landmark} value of \code{\link{img2landmark}}.}

  \item{des}{a \code{\link{routine1}} object}

  \item{nef}{a \code{\link{rNEF}} object}

  \item{gpa}{a \code{\link{rGPA}} object}

  \item{pc}{numeric. number of optimized PC range
  determined by \code{\link{pccv}}}

  \item{har}{numeric. number of optimized harmonics range
  determined by \code{\link{harcv}}}
  
  \item{class}{character/factor giving specimen's species/groups.}
}
\value{
a \code{.rds} file saved or updated.
}
\description{
save/ update the processed semi-landmarks, shape
descriptor data, GPA/ EFA transformed data, optimized
parameters into a \code{.rds} file.
}
\details{
For \code{updateproj}, the objects assigned with new values will be overwrited
with new ones, the rest will remain the same. The current version only support overwrite(!) the whole selected object, not append.

}
\seealso{
To read the project, use: \code{\link{readRDS}}

Which this function wraps: \code{\link{saveRDS}}
}