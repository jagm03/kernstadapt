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
#' @examples
#' data(lGCpp)
#' stIntensity <- dens.par(lGCpp, dimt = 16)
#' plot(spatstat.geom::as.solist(stIntensity[13:16]), ncols = 4,
#'      main = 'Non-separable Example', equal.ribbon = TRUE)
#'
#' @importFrom sparr OS
#' @importFrom stats quantile
#' @importFrom spatstat.geom setmarks marks
#' @importFrom spatstat.random rpoispp
#' @export
dens.par <- function(X, #point pattern
                     ...,
                     bw.xy = NULL, #bandwidths
                     ngroups.xy = NULL, #groups
                     at = c("pixels", "points") #at
){
  stopifnot(is.lpp(X))
  n <- npoints(X)
  at <- match.arg(at)
  ngroups.xy <- check.ngroups(ngroups.xy, N = n, order = 2)

  if (missing(bw.xy) || is.null(bw.xy)) {
    bw.xy <- bw.abram.net(X, h0 = OS(X))
  }
  else if (is.numeric(bw.xy)) {
    check.nvector(bw.xy, n, oneok = TRUE)
    if (length(bw.xy) == 1)
      bw.xy <- rep(bw.xy, n)
  }
  else stop("Argument 'bw.xy' should be a single value or a numeric vector")

  p.xy <- seq(0, 1, length = ngroups.xy + 1)
  p.t <- seq(0, 1, length = ngroups.t + 1)
  qbands.xy <- quantile(bw.xy, p.xy)
  qbands.t <- quantile(bw.t, p.t)
  groupid.xy <- findInterval(bw.xy, qbands.xy, all.inside = T)
  groupid.t <- findInterval(bw.t, qbands.t, all.inside = T)
  pmid.xy <- (p.xy[ -1 ] + p.xy[ -length(p.xy) ]) / 2
  pmid.t <- (p.t[ -1 ] + p.t[ -length(p.t) ]) / 2
  qmid.xy <- quantile(bw.xy, pmid.xy)
  qmid.t <- quantile(bw.t, pmid.t)
  group.xy <- factor(groupid.xy, levels = 1:ngroups.xy)
  group.t <- factor(groupid.t, levels = 1:ngroups.t)
  BX <- setmarks(X, data.frame(t = t, g.t = group.t))
  Y.xy <- split(BX, group.xy)

  # remove zero pp
  ok1 <- sapply(Y.xy, npoints) > 0
  Y.xy <- Y.xy[ok1]
  qmid.xy <- qmid.xy[ok1]

  Y.xy <- mapply(function(x, e) {
    u <- cbind(marks(x), epsilon = e)
    return(setmarks(x, u))
  },
  x = Y.xy, e = as.list(qmid.xy), SIMPLIFY = F)

  Y.xy.t <- unlist(lapply(Y.xy, split, f = "g.t"), recursive = F)
  Y.xy.t <- lapply(Y.xy.t, function(x) {
    x$marks$delta = qmid.t[x$marks$g.t]
    return(x)})

  PP <- lapply(Y.xy.t, function(x) setmarks(unmark(x), x$marks$t))
  B.XY <- lapply(Y.xy.t, function(x) marks(x)$epsilon[1])
  B.t <- lapply(Y.xy.t, function(x) marks(x)$delta[1])

  ok <- sapply(PP, npoints) > 0
  PP <- PP[ok]; B.XY <- B.XY[ok]; B.t <- B.t[ok]

  Z <- mapply(dens.direct.engine, X = PP, bw.xy = B.XY, bw.t = B.t, SIMPLIFY = F,
              MoreArgs = list(tlim = range.t, sres = dimyx, tres = dimt))

  ZZ <- Reduce("+", Z)
  WM <- as.mask(X$window, dimyx = dimyx)
  inside <- WM$m
  grx <- WM$xcol
  gry <- WM$yrow
  zz <- list()

  for(i in 1:dimt){
    zz[[i]] <- im(ZZ[,, i], xcol = grx, yrow = gry)
    zz[[i]][!inside] <- NA
  }

  if (at == "bins")
    return(zz)

  if (at == "points"){
    lambda <- rep(NA, nT)
    tcw <- diff(range.t) / dimt
    grt <- range.t[1] + 0.5 * tcw + (0:(dimt - 1)) * tcw
    tC <- findInterval(t, grt, all.inside = T)
    for (i in 1:nT) {
      lambda[i] <- safelookup(zz[[tC[i]]], X[i])
    }
    return(lambda)
  }
}
