\name{plotallo}
\alias{plotallo}
\title{Plot allometric deformation}
\usage{
plotallo(gpa, group, target, col = 1, main, n = 40, saveplot = FALSE,
  plotsize = 1300)
}
\arguments{
  \item{gpa}{a \code{\link{rGPA}} object}

  \item{group}{factor giving the groupings, must be in same
  order as \code{A} of \code{\link{rGPA}} object}

  \item{target}{character. group that you want to draw, one
  of the level of \code{group}}

  \item{main}{character. title of the plot. default is
  using the \code{target}}

  \item{n}{numeric. no. of grids used in TPS visualization,
  argument passed to \code{tps}}

  \item{col}{color of the outline}

  \item{saveplot}{logical. Whether to save the plot. No
  plot will be displayed in window if \code{saveplot=TRUE}}

  \item{plotsize}{width of the plot in pixel. Used only
  when \code{saveplot=TRUE}}
}
\value{
just the plot
}
\description{
A wrapper function to draw min, mean, and max shape of a
species out of all configurations of selected group. Min
and max shapes are displayed as TPS deformations from the
mean shape
}

