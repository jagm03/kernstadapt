#' Non-separable adaptive spatio-temporal intensity estimator
#'
#' Provides an adaptive-bandwidth kernel estimate for spatio-temporal point patterns in a non-separable fashion by using binning of the bandwidth values.
#'
#' @param X A spatial point pattern (an object of class \code{ppp}) with the spatial coordinates of the observations.
#' @param t A numeric vector of temporal coordinates with equal length to the number of points in \code{X}. This gives the time associated with each spatial point.
#' @param dimyx Spatial pixel resolution. The default is 128 for each axes.
#' @param dimt Temporal bin vector dimension. The default is 128.
#' @param bw.xy Numeric vector of spatial smoothing bandwidths for each point in \code{X}. By default this is computed using \link[spatstat.core]{bw.abram}.
#' @param bw.t Numeric vector of temporal smoothing bandwidths for each point in \code{t}. By default this is computed using \link{bw.abram.temp}.
#' @param ngroups.xy Number of groups in which the spatial bandwidths should be partitioned. If this number is 1, then a classical non-adaptive estimator will be used for the spatial part with a bandwidth selected as the median of the bw.xy vector.
#' @param ngroups.t Number of groups in which the temporal bandwidths should be partitioned. If this number is 1, then a classical non-adaptive estimator will be used for the temporal part with a bandwidth selected as the median of the bw.t vector.
#' @param at String specifying whether to estimate the intensity at a mesh (\code{at = "bins"}) or only at the points of \code{X} (\code{at = "points"}).
#' @details
#' This function computes a non-separable spatio-temporal adaptive kernel estimate of the intensity. It starts from a planar point pattern \code{X} and a vector of times \code{t} and partition (cells) the spatial and temporal components to apply a non-separable kernel estimator within each cell.
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
#' stIntensity <- dens.par(X,t, at = "bins")
#' plot(as.imlist(stIntensity[13:16]), main = 'Non-separable Example')
#' }
#'
#' @importFrom sparr OS
#' @importFrom stats quantile
#' @importFrom spatstat.geom setmarks marks
#' @importFrom spatstat.random rpoispp
#' @export
dens.par <- function(X, t = NULL, #point patterns
                     dimyx = 128, dimt = 128, #resolution
                     bw.xy = NULL, bw.t = NULL, #bandwidths
                     ngroups.xy = NULL, ngroups.t = NULL, #groups
                     at = c("bins", "points") #at
){
  verifyclass(X, "ppp")
  n <- npoints(X)
  if(is.null(t)) t <- marks(X)
  t <- checkt(t)
  nT <- length(t)
  if(nT != n)
    stop(paste("Length of temporal vector does not match number of spatial observations\n   npoints(X) = ",n,"; length(t) = ",length(t), sep = ""))

  at <- match.arg(at)
  range.t <- range(t)
  ngroups.xy <- check.ngroups(ngroups.xy, N = nT, order = 3)
  ngroups.t <- check.ngroups(ngroups.t, N = nT, order = 6)

  if (missing(bw.t) || is.null(bw.t)) {
    bw.t <- bw.abram.temp(t)
  }
  else if (is.numeric(bw.t)) {
    check.nvector(bw.t, nT, oneok = TRUE)
    if (length(bw.t) == 1)
      bw.t <- rep(bw.t, nT)
  }
  else stop("Argument 'bw.t' should be a single value or a numeric vector")

  if (missing(bw.xy) || is.null(bw.xy)) {
    bw.xy <- bw.abram(X, h0 = OS(X))
  }
  else if (is.numeric(bw.xy)) {
    check.nvector(bw.xy, nT, oneok = TRUE)
    if (length(bw.xy) == 1)
      bw.xy <- rep(bw.xy, nT)
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
