#' Separable adaptive spatio-temporal intensity estimator
#'
#' Provides an adaptive-bandwidth kernel estimate for spatio-temporal point patterns in a separable fashion, i.e., by multiplying spatial and temporal marginals.
#'
#' @param X A spatial point pattern (an object of class \code{ppp}) with the spatial coordinates of the observations.
#' @param t A numeric vector of temporal coordinates with equal length to the number of points in \code{X}. This gives the time associated with each spatial point.
#' @param bw.xy Numeric vector of spatial smoothing bandwidths for each point in \code{X}. By default this is computed usign \link[spatstat.core]{bw.abram}.
#' @param bw.t Numeric vector of temporal smoothing bandwidths for each point in \code{t}. By default this is computed usign \link{bw.abram.temp}.
#' @param ngroups.xy Number of groups in which the spatial bandwidths should be partitioned. If this number is 1, then a classical non-adaptive estimator will be used for the spatial part with a bandwidth selected as the median of the bw.xy vector.
#' @param ngroups.t Number of groups in which the temporal bandwidths should be partitioned. If this number is 1, then a classical non-adaptive estimator will be used for the temporal part with a bandwidth selected as the median of the bw.t vector.
#' @param dimt Temporal bin vector dimension. The default is 128.
#' @param dimyx Spatial pixel resolution. The default is 128 for each axes.
#' @param at String specifying whether to estimate the intensity at a mesh (\code{at = "bins"}) or only at the points of \code{X} (\code{at = "points"}).
#'
#' @details
#' This function computes a spatio-temporal adaptive kernel estimate of the intensity in a separable fashion. It starts from a planar point pattern \code{X} and a vector of times \code{t} and uses partition techniques separately for the spatial and temporal components, then it multiplies both components and normalises by the number of points to preserve the mass.
#' The arguments \code{bw.xy} and \code{bw.t} specify the smoothing bandwidth vectors to be applied to each of the points in \code{X} and \code{t}. They should be a numeric vectors of bandwidths.
#' The method partition the range of bandwidths into intervals, subdividing the points of the pattern \code{X} and \code{t} into sub-patterns according to the bandwidths, and applying fixed-bandwidth smoothing to each sub-pattern. Specifying \code{ngroups.xy = 1} is the same as fixed-bandwidth smoothing with bandwidth \code{sigma = median(bw.xy)} in the spatial case and \code{ngroups.t = 1} is the same as fixed-bandwidth smoothing with bandwidth \code{sigma = median(bw.xy)}.
#'
#' @return
#' If \code{at = "points"} (the default), the result is a numeric vector with one entry for each data point in \code{X}. if \code{at = "bins"} is a list named (by time-point) list of pixel images (\link[spatstat.geom]{im} objects) corresponding to the joint spatio-temporal intensity over space at each discretised time bin.
#'
#' @author Jonatan A. González
#'
#' @examples
#' \dontrun{
#' X <- rpoispp(1400)
#' t <- rbeta(X$n, 1,4,0.8)
#' stIntensity <- dens.st.sep.bin(X,t, at = "bins")
#' plot(as.imlist(stIntensity[13:16]), main = 'Separable Example')
#' }
#'
#' @importFrom spatstat.core densityAdaptiveKernel.ppp
#' @importFrom spatstat.geom eval.im
#' @export
dens.st.sep.bin <- function(X, t = NULL,
                            bw.xy = NULL, bw.t = NULL, #bandwidths
                            ngroups.xy = NULL, ngroups.t = NULL, #groups
                            dimt = 128, dimyx = 128, #resolution
                            at = c("bins", "points") #at
)
{
  verifyclass(X, "ppp")
  n <- npoints(X)
  if(is.null(t)) t <- marks(X)
  t <- checkt(t)
  if(length(t) != n)
    stop(paste("Length of temporal vector does not match number of spatial observations\n   npoints(X) = ",n,"; length(t) = ",length(t), sep = ""))

  at <- match.arg(at)
  at.s <- switch (at, bins = "pixels", points = "points")
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

  Mtemporal <- dens.t.bin(ti = t, bw = bw.t,
                          ngroups = ngroups.t,
                          dimt = dimt, at = at)
  Mspatial <- densityAdaptiveKernel.ppp(unmark(X), bw = bw.xy, at = at.s,
                                        dimyx = dimyx, edge = T,
                                        ngroups = ngroups.xy)
  Mst <- switch(at,
                points = Mtemporal * Mspatial / length(t),
                bins = lapply(Mtemporal$y, function(x) eval.im(Mspatial * x / length(t)))
  )
  return(Mst)
}
