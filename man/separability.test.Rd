\name{separability.test}
\alias{separability.test}
\title{
Separability test for spatio-temporal point patterns
}
\description{
Performs a separability test of the first-order intensity function based on Fisher's test and counts.
}
\usage{
separability.test(X, t,
                  nx = NULL, ny = NULL, nt = NULL,
                  nperm = 1000)
}
\arguments{
  \item{X}{A spatial point pattern (an object of class \code{ppp}) with the spatial coordinates of the observations}
  \item{t}{A numeric vector of temporal coordinates with equal length to the number of points in \code{X}. This gives the time associated with each spatial point.}
 \item{nx, ny, nt}{Numbers of quadrats in the \eqn{x,y} and \eqn{t} directions.}
 \item{nperm}{An integer specifying the number of replicates used in the Monte Carlo test.}
}
\details{
This function computes a non-separable spatio-temporal adaptive kernel estimate of the intensity. It starts from a planar point pattern \code{X} and a vector of times \code{t} and partition (cells) the spatial and temporal components to apply a non-separable kernel estimator within each cell.
The arguments \code{bw.xy} and \code{bw.t} specify the smoothing bandwidth vectors to be applied to each of the points in \code{X} and \code{t}. They should be a numeric vectors of bandwidths.

The method partition the range of bandwidths into intervals, subdividing the points of the pattern \code{X} and \code{t} into sub-patterns according to the bandwidths, and applying fixed-bandwidth smoothing to each sub-pattern. Specifying \code{ngroups.xy = 1} is the same as fixed-bandwidth smoothing with bandwidth \code{sigma = median(bw.xy)} in the spatial case and \code{ngroups.t = 1} is the same as fixed-bandwidth smoothing with bandwidth \code{sigma = median(bw.xy)}.
}
\value{
If \code{at = "points"} (the default), the result is a numeric vector with one entry for each data point in \code{X}. if \code{at = "bins"} is a list named (by time-point) list of pixel images (\link[spatstat.geom]{im} objects) corresponding to the joint spatio-temporal intensity over space at each discretised time bin.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Jonatan A. González
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
X <- rpoispp(1400)
t <- rbeta(X$n, 1,4,0.8)
stIntensity <- dens.st.bin(X,t, at = "bins")
plot(as.imlist(stIntensity[13:16]), main = 'Non-separable Example')
}
}
