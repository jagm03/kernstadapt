#' Direct separable adaptive spatio-temporal intensity estimator
#'
#' Provides an adaptive-bandwidth kernel estimate for spatio-temporal point patterns in a separable fashion, i.e., by multiplying spatial and temporal marginals. This estimation is performed by calculating the classical estimator, i.e., the slowest estimation.
#'
#' @param X A spatial point pattern (an object of class \code{ppp}) with the spatial coordinates of the observations. It may contain marks representing times.
#' @param t A numeric vector of temporal coordinates with equal length to the number of points in \code{X}. This gives the time associated with each spatial point. This argument is not necessary if time marks are provided to the point pattern \code{X}.
#' @param dimyx Spatial pixel resolution. The default is 128 for each axes.
#' @param dimt Temporal bin vector dimension. The default is 128.
#' @param bw.xy Numeric vector of spatial smoothing bandwidths for each point in \code{X}. By default this is computed using \link[spatstat.core]{bw.abram}.
#' @param bw.t Numeric vector of temporal smoothing bandwidths for each point in \code{t}. By default this is computed using \link{bw.abram.temp}.
#' @param at String specifying whether to estimate the intensity at a mesh (\code{at = "bins"}) or only at the points of \code{X} (\code{at = "points"}).
#'
#' @details
#' This function computes a spatio-temporal adaptive kernel estimate of the intensity in a separable fashion. It starts from a planar point pattern \code{X} and a vector of times \code{t} and uses a direct estimator for each dimension, then it multiplies both components and normalises by the number of points to preserve the mass.
#' The arguments \code{bw.xy} and \code{bw.t} specify the smoothing bandwidth vectors to be applied to each of the points in \code{X} and \code{t}. They should be a numeric vectors of bandwidths.
#' The method partition the range of bandwidths into intervals, subdividing the points of the pattern \code{X} and \code{t} into sub-patterns according to the bandwidths, and applying fixed-bandwidth smoothing to each sub-pattern. Specifying \code{ngroups.xy = 1} is the same as fixed-bandwidth smoothing with bandwidth \code{sigma = median(bw.xy)} in the spatial case and \code{ngroups.t = 1} is the same as fixed-bandwidth smoothing with bandwidth \code{sigma = median(bw.xy)}.
#'
#' @return
#' If \code{at = "points"}, the result is a numeric vector with one entry for each data point in \code{X}. if \code{at = "bins"} is a list named (by time-point) list of pixel images (\link[spatstat.geom]{im} objects) corresponding to the joint spatio-temporal intensity over space at each discretised time bin.
#'
#' @references
#' González J.A. and Moraga P. (2018)
#' An adaptive kernel estimator for the intensity function of spatio-temporal point processes
#' <https://arxiv.org/pdf/2208.12026.pdf>
#'
#' @author Jonatan A. González
#'
#' @examples
#' data(santander)
#' X <- santander[sample.int(200)]
#' stIntensity <- dens.direct.sep(X,
#'                                dimyx = 32, dimt = 6,
#'                                at = "bins")
#' plot(spatstat.geom::as.solist(stIntensity[2:4]),
#'      main = 'Direct separable example', equal.ribbon = TRUE)
#'
#' @export
dens.direct.sep <- function(X, t = NULL,
                            dimyx = 128, dimt = 128, #resolution
                            bw.xy = NULL, bw.t = NULL, #bandwidths
                            at = c("bins", "points") #at
)
{
  verifyclass(X, "ppp")
  n <- npoints(X)
  if(is.null(t)) t <- marks(X)
  t <- checkt(t)
  nT <- length(t)
  if(length(t) != n)
    stop(paste("Length of temporal vector does not match number of spatial observations\n   npoints(X) = ",n,"; length(t) = ",length(t), sep = ""))
  at <- match.arg(at)
  at.s <- switch (at, bins = "pixels", points = "points")
  range.t <- range(t)
  if (is.null(bw.xy)) {
    h0 <- OS(X)
    bw.xy <- bw.abram(X, h0)
  }
  if (missing(bw.t) || is.null(bw.t)) {
    bw.t <- bw.abram.temp(t)
  }
  else if (is.numeric(bw.t)) {
    check.nvector(bw.t, nT, oneok = TRUE)
    if (length(bw.t) == 1)
      bw.t <- rep(bw.t, nT)
  }
  #Temporal part
  PPt <- split(t, f = factor(1:X$n))
  Z <- mapply(density.default, x = PPt, bw = bw.t, SIMPLIFY = F,
              MoreArgs = list(kernel = "gaussian", n = dimt,
                              from = range.t[1], to = range.t[2] ))
  zx <- Z[[1]]$x
  Edge0 <- sapply(bw.t, edge.t, dt = zx, TS = range.t[2], TI = range.t[1])
  Ys <- sapply(Z, "[[", 'y') / Edge0
  Ys <- apply(Ys, 1, sum)
  Mtemporal <- switch(at,
                      points = approxfun(x = zx, y = Ys)(t),
                      bins = Ys)

  #Spatial part
  at.s <- switch (at, bins = "pixels", points = "points")
  Mspatial <- densityAdaptiveKernel.ppp(unmark(X), bw = bw.xy, at = at.s,
                                        dimyx = dimyx, edge = T, ngroups = Inf)
  Mst <- switch(at,
                points = Mtemporal * Mspatial / length(t),
                bins = lapply(Mtemporal, function(x) eval.im(Mspatial * x / length(t)))
  )
  return(Mst)
}
