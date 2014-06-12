\name{otolith-package}
\alias{otolith-package}
\alias{otolith}
\docType{package}
\title{
Identification of fish otolith using digital images
}
\description{
A package provides tools to process digital images of fish otoliths and building identification tools based on geometric morphometrics.

This package forms a part of the project on the study of Malaysian fish otolith with the following members:
  \itemize{
    \item Wong Jin Yung \email{wongjinyung@gmail.com}
    \item Loh Kar Hoe \email{khloh@um.edu.my}
    \item Cecilia Chu \email{coccynx@yahoo.com}
    \item Chong Ving Ching \email{chong@um.edu.my}
  }
}
\details{
\emph{Featured functions in this package:}
  
  \strong{Image processing}
  \tabular{ll}{
    \code{\link{extractout}} \tab automated outline extraction from otolith image 
    \cr 
    \code{\link{comb}} / \code{\link{equaldist}} \tab automated sampling of 
    curvilinear semi-landmarks \cr
    \code{\link{img2landmark}} \tab wrapper to batch convert otolith images into 
    semi-landmarks \cr
  }
  \strong{Classification evaluation}
  \tabular{ll}{
    \code{\link{mrkfcv}} / \code{\link{mrkfcv2}} \tab iterated k-fold cross
    -validation (kfcv)\cr
    \code{\link{pccv}} / \code{\link{harcv}} / \code{\link{threcv}} \tab
    evaluation of dimension reduction / posterior probability thresholding based 
    on evalution of overall accuracy through iterated kfcv \cr
  }
  \strong{Search and Prediction}
  \tabular{ll}{
    \code{\link{otosearch}} \tab search unknown otolith image(s) against a project 
    database, regardless of otolith's side and orientation \cr
    \code{\link{otopred}} \tab predict unknown otolith image(s) against a project 
    database, based on classification model \cr
  }
}