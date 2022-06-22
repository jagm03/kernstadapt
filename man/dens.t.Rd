\name{dens.t.bin}
\alias{dens.t.bin}
\title{
Adaptive kernel estimate of intensity on a temporal point pattern
}
\description{
Estimates the intensity of a point process with only temporal dimension by applying
an adaptive (variable bandwidth) Gaussian edge-corrected kernel smoothing.
}
\usage{
dens.t.bin(ti, bw = NULL, ngroups = NULL, dimt = 128, at = c("points", "bins"))
}
\arguments{
  \item{ti}{Temporal point pattern, a vector with observations.}
  \item{bw}{Numeric vector of smoothing bandwidths for each point in ti. The default is to compute bandwidths using \link{bw.abram.temp}.}
  \item{ngroups}{Number of groups into which the bandwidths should be partitioned and discretised. The default is the square root (rounded) of the number of points of \code{ti}.}
  \item{dimt}{Bin vector dimension. The default is 128.}
  \item{at}{String specifying whether to estimate the intensity at bins points (\code{at = "bins"}) or only at the points of ti (\code{at = "points"}).}
}
\details{
This function computes a temporally-adaptive kernel estimate of the intensity from a one-dimensional point pattern ti using the partitioning technique of Davies and Baddeley (2018).
The argument bw specifies the smoothing bandwidths to be applied to each of the points in X. It should be a numeric vector of bandwidths.
Let the points of \eqn{ti} be \eqn{t_1, ..., t_n} and the corresponding bandwidths \eqn{\sigma_1,...,\sigma_n}, then the adaptive kernel estimate of intensity at a location \eqn{t} is

\deqn{\lambda(t) = \sum_{i=1}^n \frac{K(t,t_i; \sigma_i)}{c(t; \sigma_i)}}

where \eqn{K()} is the Gaussian smoothing kernel.

The method partition the range of bandwidths into ngroups intervals, correspondingly subdividing the points of the pattern \code{ti} into ngroups sub-patterns according to bandwidth, and applying fixed-bandwidth smoothing to each sub-pattern. Specifying \code{ngroups = 1} is the same as fixed-bandwidth smoothing with bandwidth \code{sigma = median(bw)}.
}
\value{
If \code{at = "points"} (the default), the result is a numeric vector with one entry for each data point in \code{ti}. If \code{at = "bins"} the result is a data.frame containing the \eqn{x,y} coordinates of the intensity function.
}
\references{
Davies, T.M. and Baddeley, A. (2018) Fast computation of spatially adaptive kernel estimates. Statistics and Computing, 28(4), 937-956.
}
\author{
Jonatan A. González
}
\examples{
ti <- rbeta(100, 1,4,0.8)
dens.t.bin(ti)
}
