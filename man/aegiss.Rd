\name{aegiss}
\alias{aegiss}
\docType{data}
\title{
Diggle et al.'s data: Non-specific gastrointestinal data
}
\description{
A spatio-temporal point pattern of locations and times of individuals with non-specific gastrointestinal infections in Hampshire, UK, from 2001 to 2003.
}
\usage{data("aegiss")}
\format{
  An object of class \code{"ppp"}
  Entries include
  \tabular{ll}{
    \code{x} \tab Cartesian \eqn{x}-coordinate of infection \cr
    \code{y} \tab Cartesian \eqn{y}-coordinate of infection \cr
    \code{marks} \tab An integer vector of values indicating \cr
                      the time of infection starting in \code{1}
    \code{window} \tab Cartesian coordinates of Hampshire's map
    }
  See \code{\link[spatstat.geom]{ppp.object}} for details of the format.
}
\source{
Diggle, P. (2015) AEGISS1. Syndromic surveillance of gastro-intestinal illness
\url{https://www.research.lancs.ac.uk/portal/en/datasets/aegiss1-syndromic-surveillance-of-gastrointestinal-illness(e732159d-811f-44b0-b2d3-0139ebf4edee).html}
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
