#' Non-separable adaptive spatio-temporal intensity estimator
#'
#' Provides an adaptive-bandwidth kernel estimate for spatio-temporal point patterns in a non-separable fashion by using binning of the bandwidth values.
#'
#' @param X A point pattern on a linear network (an object of class \code{lpp}) to be smoothed.
#' @param bw.xy Numeric vector of spatial smoothing bandwidths for each point in \code{X}. By default this is computed using \link[spatstat.explore]{bw.abram}.
#' @param ngroups.xy Number of groups in which the bandwidths should be partitioned. If this number is 1, then a classical non-adaptive estimator will be used for the spatial part with a bandwidth selected as the median of the bw.xy vector.
#' @param at String specifying whether to estimate the intensity at a mesh (\code{at = "pixels"}) or only at the points of \code{X} (\code{at = "points"}).
#' @details
#' This function computes an adaptive kernel estimate of the intensity on linear networks. It starts from a point pattern \code{X} and partition the spatial component to apply a kernel estimator within each cell.
#' The argument \code{bw.xy} specify the smoothing bandwidth vectors to be applied to each of the points in \code{X}. It should be a numeric vector of bandwidths.
#' The method partition the range of bandwidths into intervals, subdividing the points of the pattern \code{X} into sub-patterns according to the bandwidths, and applying fixed-bandwidth smoothing to each sub-pattern. Specifying \code{ngroups.xy = 1} is the same as fixed-bandwidth smoothing with bandwidth \code{sigma = median(bw.xy)} in the spatial case and \code{ngroups.t = 1} is the same as fixed-bandwidth smoothing with bandwidth \code{sigma = median(bw.xy)}.
#'
#' @return
#' If \code{at = "points"} (the default), the result is a numeric vector with one entry for each data point in \code{X}. if \code{at = "pixels"} is a pixel image on a linear network (\link[spatstat.linnet]{linim} objects) corresponding to the intensity over linear network.
#'
#' @references
#' González J.A. and Moraga P. (2018)
#' An adaptive kernel estimator for the intensity function of spatio-temporal point processes
#' <https://arxiv.org/pdf/2208.12026.pdf>
#'
#' @author Jonatan A. González
#'
#'
#' @importFrom sparr OS
#' @importFrom stats quantile
#' @importFrom spatstat.geom setmarks marks
#' @importFrom spatstat.random rpoispp
#' @export
dens.net.heat <- function(X, #point pattern
                          ...,
                          weights = NULL, #optional weights
                          bw = NULL, #bandwidths
                          ngroups = NULL, #groups
                          at = c("pixels", "points") #at
){
  stopifnot(is.lpp(X))
  at <- match.arg(at)
  nX <- npoints(X)
  ngroups <- check.ngroups(ngroups, N = nX, order = 2)

  if(weighted <- !is.null(weights)) {
    check.nvector(weights, nX, oneok = TRUE, vname = "weights")
    if(length(weights) == 1) weights <- rep(weights, nX)
  } else weights <- rep(1, nX)

  if (missing(bw) || is.null(bw)) {
    bw <- bw.abram.net(X, h0 = OS(X))
  }
  else if (is.numeric(bw)) {
    check.nvector(bw, nX, oneok = TRUE)
    if (length(bw) == 1)
      bw <- rep(bw, nX)
  }
  else stop("Argument 'bw' should be a single value or a numeric vector")

  #divide bandwidths into groups
  if(ngroups == nX) {
    ## every data point is a separate group
    groupid <- 1:nX
    qmid <- bw
  } else {
    # usual case
    p <- seq(0, 1, length = ngroups + 1)
    qbands <- quantile(bw, p)
    groupid <- findInterval(bw, qbands, all.inside = TRUE)
    # map to middle of group
    pmid <- (p[-1] + p[-length(p)]) / 2
    qmid   <- quantile(bw, pmid)
  }

  marks(X) <- if(weighted) weights else NULL
  group <- factor(groupid, levels = 1:ngroups)
  Y <- split(X, group)

  Z <- mapply(densityHeat.lpp,
              x = Y,
              sigma = as.list(qmid),
              SIMPLIFY = F,
              MoreArgs = list(at = at, ...))

  ZZ <- switch(at,
               pixels = im.apply(Z, "sum"),
               points = unsplit(Z, group.xy))
  return(ZZ)
}
