\name{santander}
\alias{santander}
\docType{data}
\title{
Santander earthquakes epicentres
}
\description{
A spatio-temporal point pattern of locations and times of earthquake epicentres in Santander-Colombia with a magnitude above 3.5ML, recorded from January 2010 to December 2020.
}
\usage{data("santander")}
\format{
  An object of class \code{"ppp"}
  Entries include
  \tabular{ll}{
    \code{x} \tab Cartesian \eqn{x}-coordinate of infection \cr
    \code{y} \tab Cartesian \eqn{y}-coordinate of infection \cr
    \code{marks} \tab An integer vector of values indicating \cr
                      the time of the earthquake starting in \code{1}
    \code{window} \tab Cartesian cordinates of Santander's map
    }
  See \code{\link[spatstat.geom]{ppp.object}} for details of the format.
}
\source{
González et. al., (2021) Classification of Events Using Local Pair Correlation Functions for Spatial Point Patterns. \emph{Journal of Agricultural, Biological and Environmental Statistics.} 26(4) \bold{538-559.}
}
\references{
Diggle, L. Knorr-Held, R. B, T. Su, P. Hawtin, and T. Bryant. (2003)
On-line monitoring of public health surveillance data.
In B. R. and S. D.F., editors, Monitoring the Health of Populations:
\emph{Statistical Principles and Methods for Public Health Surveillance.} \bold{233–266}.
Oxford University Press.
}
\examples{
data(aegiss)
head(aegiss)
plot(aegiss, bg = terrain.colors(20), pch = 21, cex = 1)
}
\keyword{datasets}
\keyword{aegiss}
