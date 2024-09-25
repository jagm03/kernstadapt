#' Adaptive linear netwokr intensity estimator using the Okabe-Sugihara equal-split algorithms
#'
#' Computes an adaptive-bandwidth kernel estimate for the intensity function through the Okabe-Sugihara equal-split algorithms by using binning of the bandwidth values.
#'
#' @param X A point pattern on a linear network (an object of class \code{lpp}) to be smoothed.
#' @param ... Extra arguments passed to \link[spatstat.linnet]{densityHeat.lpp}.
#' @param weights Optional. Numeric vector of weights associated with the points of X. Weights can be positive, negative or zero.
#' @param bw Numeric vector of spatial smoothing bandwidths for each point in \code{X}. By default this is computed using \link[spatstat.univar]{bw.abram}.
#' @param ngroups Number of groups in which the bandwidths should be partitioned. If this number is 1, then a classical non-adaptive estimator will be used for the spatial part with a bandwidth selected as the median of the bw.xy vector.
#' @param at String specifying whether to estimate the intensity at a mesh (\code{at = "pixels"}) or only at the points of \code{X} (\code{at = "points"}).
#' @param verbose Logical value indicating whether to print progress reports for every partition group.
#' @details
#' This function computes an adaptive kernel estimate of the intensity on linear networks. It starts from a point pattern \code{X} and partition the spatial component to apply a kernel estimator within each cell.
#' The argument \code{bw} specifies the smoothing bandwidth vector to be applied to each of the points in \code{X}. It should be a numeric vector of bandwidths.
#' The method partition the range of bandwidths into intervals, subdividing the points of the pattern \code{X} into sub-patterns according to the bandwidths, and applying fixed-bandwidth smoothing to each sub-pattern. Specifying \code{ngroups = 1} is the same as fixed-bandwidth smoothing with bandwidth \code{sigma = median(bw)}.
#'
#' @return
#' If \code{at = "points"} (the default), the result is a numeric vector with one entry for each data point in \code{X}. if \code{at = "pixels"} is a pixel image on a linear network (\link[spatstat.linnet]{linim} objects) corresponding to the intensity over linear network.
#'
#' @references
#' González J.A. and Moraga P. (2018)
#' An adaptive kernel estimator for the intensity function of spatio-temporal point processes
#' <http://arxiv.org/pdf/2208.12026>
#'
#' @author Jonatan A. González
#'
#'
#' @importFrom spatstat.linnet densityEqualSplit
#' @export
dens.net.EqualSplit <- function(X, #point pattern on a linear network
                                ...,
                                weights = NULL, #optional weights
                                bw = NULL, #bandwidths
                                ngroups = NULL, #groups
                                at = c("pixels", "points"), #at
                                verbose = FALSE
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
    bw <- bw.abram.net(X, h0 = OS(as.ppp(X)))
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

  X <- setmarks(X, if(weighted) weights else NULL)
  group <- factor(groupid, levels = 1:ngroups)
  Y <- split(X, group)

  Z <- mapply(densityEqualSplit,
              x = Y,
              sigma = as.list(qmid),
              SIMPLIFY = F,
              MoreArgs = list(at = at, verbose = verbose, ...))

  ZZ <- switch(at,
               pixels = as.linim(im.apply(Z, "sum"), L = X$domain),
               points = unsplit(Z, group))
  return(ZZ)
}
