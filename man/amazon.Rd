\name{amazon}
\alias{aegiss}
\docType{data}
\title{
Locations of fires in the Amazon biome.
}
\description{
A spatio-temporal point pattern of locations and times of active deforestation fires (starting within past 24 hours) from 01/01/2021 to 10/10/2021 (284 days).
}
\usage{data("amazon")}
\format{
  An object of class \code{"ppp"}
  Entries include
  \tabular{ll}{
    \code{x} \tab Cartesian \eqn{x}-coordinate of fire \cr
    \code{y} \tab Cartesian \eqn{y}-coordinate of fire \cr
    \code{marks} \tab An integer vector of values indicating \cr
                      the time of the fire starting in \code{2}
    \code{window} \tab Coordinates of Amazonia biome map
    }
  See \code{\link[spatstat.geom]{ppp.object}} for details of the format.
}
\source{
The Amazon Dashboard (\url{http://www.globalfiredata.org})
}
\examples{
\dontrun{
data(amazon)
head(amazon)
# Plot a sample
plot(amazon[sample(amazon$n, 1000)], bg = heat.colors(20), pch = 21, cex = 1)
}
}
\keyword{datasets}
\keyword{aegiss}
