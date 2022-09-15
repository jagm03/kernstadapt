\name{lGCpp}
\alias{lGCpp}
\docType{data}
\title{
Simulated spatio-temporal log Gaussian Cox point pattern data
}
\description{
A simulated spatio-temporal log Gaussian Cox point pattern with a underlying non-separable Geneiting covariance.
}
\usage{data("lGCpp")}
\format{
  An object of class \code{"ppp"}
  Entries include
  \tabular{ll}{
    \code{x} \tab Cartesian \eqn{x}-coordinates \cr
    \code{y} \tab Cartesian \eqn{y}-coordinates \cr
    \code{marks} \tab A vector of values indicating the times \cr
    \code{window} \tab Cartesian coordinates of the unit square
    }
  See \code{\link[spatstat.geom]{ppp.object}} for details of the format.
}

\references{
Gabriel E, Rowlingson B, Diggle PJ (2013). \code{stpp:} An R Package for Plotting, Simulating and Analyzing Spatio-Temporal Point Patterns. \emph{Journal of Statistical Software,} \bold{(53)(2)}(2), 1–29. \cr

González J.A. and Moraga P. (2018)
An adaptive kernel estimator for the intensity function of spatio-temporal point processes.
<arXiv:2208.12026>
}
\examples{
data(lGCpp)
head(lGCpp)
plot(lGCpp, bg = rainbow(250), pch = 21, cex = 1)
}
\keyword{datasets}
\keyword{lGCpp}
