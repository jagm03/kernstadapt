% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dens.direct.sep.R
\name{dens.direct.sep}
\alias{dens.direct.sep}
\title{Direct separable adaptive spatio-temporal intensity estimator}
\usage{
dens.direct.sep(
  X,
  t = NULL,
  dimt = 128,
  dimyx = 128,
  bw.xy = NULL,
  bw.t = NULL,
  at = c("bins", "points")
)
}
\arguments{
\item{X}{A spatial point pattern (an object of class \code{ppp}) with the spatial coordinates of the observations.}

\item{t}{A numeric vector of temporal coordinates with equal length to the number of points in \code{X}. This gives the time associated with each spatial point.}

\item{dimt}{Temporal bin vector dimension. The default is 128.}

\item{dimyx}{Spatial pixel resolution. The default is 128 for each axes.}

\item{bw.xy}{Numeric vector of spatial smoothing bandwidths for each point in \code{X}. By default this is computed usign \link[spatstat.core]{bw.abram}.}

\item{bw.t}{Numeric vector of temporal smoothing bandwidths for each point in \code{t}. By default this is computed usign \link{bw.abram.temp}.}

\item{at}{String specifying whether to estimate the intensity at a mesh (\code{at = "bins"}) or only at the points of \code{X} (\code{at = "points"}).}
}
\value{
If \code{at = "points"}, the result is a numeric vector with one entry for each data point in \code{X}. if \code{at = "bins"} is a list named (by time-point) list of pixel images (\link[spatstat.geom]{im} objects) corresponding to the joint spatio-temporal intensity over space at each discretised time bin.
}
\description{
Provides an adaptive-bandwidth kernel estimate for spatio-temporal point patterns in a separable fashion, i.e., by multiplying spatial and temporal marginals. This estimation is performed by calculating the classical estimator, i.e., the slowest estimation.
}
\details{
This function computes a spatio-temporal adaptive kernel estimate of the intensity in a separable fashion. It starts from a planar point pattern \code{X} and a vector of times \code{t} and uses a direct estimator for each dimension, then it multiplies both components and normalises by the number of points to preserve the mass.
The arguments \code{bw.xy} and \code{bw.t} specify the smoothing bandwidth vectors to be applied to each of the points in \code{X} and \code{t}. They should be a numeric vectors of bandwidths.
The method partition the range of bandwidths into intervals, subdividing the points of the pattern \code{X} and \code{t} into sub-patterns according to the bandwidths, and applying fixed-bandwidth smoothing to each sub-pattern. Specifying \code{ngroups.xy = 1} is the same as fixed-bandwidth smoothing with bandwidth \code{sigma = median(bw.xy)} in the spatial case and \code{ngroups.t = 1} is the same as fixed-bandwidth smoothing with bandwidth \code{sigma = median(bw.xy)}.
}
\examples{
\dontrun{
X <- rpoispp(1400)
t <- rbeta(X$n, 1,4,0.8)
stIntensity <- dens.direct.sep(X, t, at = "bins")
plot(as.imlist(stIntensity[13:16]), main = 'Direct separable example')
}

}
\author{
Jonatan A. González
}