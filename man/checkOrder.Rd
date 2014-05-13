\name{checkOrder}
\alias{checkOrder}
\title{Check outline order}
\usage{
checkOrder(test)
}
\arguments{
  \item{test}{a matrix consists of x and y coordinates, the
  outline data to be checked.}
}
\value{
logical. Whether the outline is ordered.
}
\description{
Checking whether the outline XY coordinates are in order,
i.e.  whether successive points' indices are in successive
order?  This is important because 'equaldist' function
assumes the coordinates are ordered.
}

