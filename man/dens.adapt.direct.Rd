\name{dens.adapt.direct}
\alias{dens.adapt.direct}
\title{
Direct estimation of non-separable adaptive spatio-temporal intensity estimator
}
\description{
Provides an adaptive-bandwidth kernel estimate for spatio-temporal point patterns in a non-separable fashion by calculating the classical estimator, i.e., the slowest estimation.
}
\usage{
dens.adapt.direct(X, t,
                  bw.xy = NULL, bw.t = NULL,
                  dimt = 128, dimyx = 128,
                  at = c("points", "bins"))
}
\arguments{
  \item{X}{A spatial point pattern (an object of class \code{ppp}) with the spatial coordinates of the observations}
  \item{t}{A numeric vector of temporal coordinates with equal length to the number of points in \code{X}. This gives the time associated with each spatial point.}
  \item{bw.xy}{Numeric vector of spatial smoothing bandwidths for each point in \code{X}. By default this is computed usign \link[spatstat.core]{bw.abram}.}
  \item{bw.t}{Numeric vector of temporal smoothing bandwidths for each point in \code{t}. By default this is computed usign \link{bw.abram.temp}.}
  \item{dimt}{Temporal bin vector dimension. The default is 128.}
  \item{dimyx}{Spatial pixel resolution. The default is 128 for each axes.}
  \item{at}{String specifying whether to estimate the intensity at a mesh (\code{at = "bins"}) or only at the points of \code{X} (\code{at = "points"}).}
}
\details{
This function computes a non-separable spatio-temporal adaptive kernel estimate of the intensity. It starts from a planar point pattern \code{X} and a vector of times \code{t} and apply a non-separable kernel estimator for each of the points of \code{X}. The arguments \code{bw.xy} and \code{bw.t} specify the smoothing bandwidth vectors to be applied to each of the points in \code{X} and \code{t}. They should be a numeric vectors of bandwidths.
}
\value{
If \code{at = "points"} (the default), the result is a numeric vector with one entry for each data point in \code{X}. if \code{at = "bins"} is a list named (by time-point) list of pixel images (\link[spatstat.geom]{im} objects) corresponding to the joint spatio-temporal intensity over space at each discretised time bin.
}
\author{
Jonatan A. González
}
\examples{
\dontrun{
X <- rpoispp(400)
t <- rbeta(X$n, 1,4,0.8)
stIntensity <- dens.adapt.direct(X,t, at = "bins", dimt = 32)
plot(as.imlist(stIntensity[13:16]), main = 'Non-separable direct example')
}
}
