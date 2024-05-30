#' Adaptive kernel estimate of intensity of a temporal point pattern
#'
#' Estimates the intensity of a point process with only temporal dimension by applying an adaptive (variable bandwidth) Gaussian edge-corrected kernel smoothing.
#'
#' @param t Temporal point pattern, a vector with observations.
#' @param dimt Bin vector dimension. The default is 128.
#' @param bw.t Numeric vector of smoothing bandwidths for each point in t. The default is to compute bandwidths using \link{bw.abram.temp}.
#' @param ngroups.t Number of groups into which the bandwidths should be partitioned and discretised. The default is the square root (rounded) of the number of points of \code{t}.
#' @param at String specifying whether to estimate the intensity at bins points (\code{at = "bins"}) or only at the points of t (\code{at = "points"}).
#'
#' @details
#' This function computes a temporally-adaptive kernel estimate of the intensity from a one-dimensional point pattern t using the partitioning technique of Davies and Baddeley (2018).
#' The argument bw.t specifies the smoothing bandwidths to be applied to each of the points in X. It should be a numeric vector of bandwidths.
#' Let the points of \eqn{t} be \eqn{t_1, ..., t_n} and the corresponding bandwidths \eqn{\sigma_1,...,\sigma_n}, then the adaptive kernel estimate of intensity at a location \eqn{v} is
#' \deqn{\lambda(v) = \sum_{i=1}^n \frac{K(v,t_i; \sigma_i)}{c(t; \sigma_i)}}
#' where \eqn{K()} is the Gaussian smoothing kernel.
#' The method partition the range of bandwidths into ngroups.t intervals, correspondingly subdividing the points of the pattern \code{t} into ngroups.t sub-patterns according to bandwidth, and applying fixed-bandwidth smoothing to each sub-pattern. Specifying \code{ngroups.t = 1} is the same as fixed-bandwidth smoothing with bandwidth \code{sigma = median(bw.t)}.
#'
#' @return
#' If \code{at = "points"} (the default), the result is a numeric vector with one entry for each data point in \code{t}. If \code{at = "bins"} the result is a data.frame containing the \eqn{x,y} coordinates of the intensity function.
#'
#' @references
#' Davies, T.M. and Baddeley, A. (2018) Fast computation of spatially adaptive kernel estimates. Statistics and Computing, 28(4), 937-956.
#'
#' González J.A. and Moraga P. (2018)
#' An adaptive kernel estimator for the intensity function of spatio-temporal point processes
#' <http://arxiv.org/pdf/2208.12026>
#'
#' @author
#' Jonatan A. González
#'
#' @examples
#' t <- rbeta(100, 1,4,0.8)
#' tIntensity <- dens.par.temp(t, at = "bins")
#' plot(tIntensity$x, tIntensity$y, type = "l")
#'
#' @export
dens.par.temp <- function(t,
                          dimt = 128,
                          bw.t = NULL,
                          ngroups.t = NULL,
                          at = c("bins", "points")){
  stopifnot(sum(t < 0) == 0)
  at <- match.arg(at)
  nT <- length(t)
  range.t <- range(t)
  if (missing(ngroups.t) || is.null(ngroups.t)) {
    ngroups.t <- max(1L, floor(sqrt(nT)))
  }
  else if (any(is.infinite(ngroups.t))) {
    ngroups.t <- nT
  }
  else {
    check.1.integer(ngroups.t)
    ngroups.t <- min(nT, ngroups.t)
  }

  if (missing(bw.t) || is.null(bw.t)) {
    bw.t <- bw.abram.temp(t)
  }
  else if (is.numeric(bw.t)) {
    check.nvector(bw.t, nT, oneok = TRUE)
    if (length(bw.t) == 1)
      bw.t <- rep(bw.t, nT)
  }
  else stop("Argument 'bw.t' should be a single value or a numeric vector")

  p <- seq(0, 1, length = ngroups.t + 1)
  qbands <- quantile(bw.t, p)
  groupid <- findInterval(bw.t, qbands, all.inside = T)
  pmid <- (p[-1] + p[ -length(p) ]) / 2
  qmid <- quantile(bw.t, pmid)
  group <- factor(groupid, levels = 1:ngroups.t)
  Y <- split(t, group, drop = F)
  ni <- sapply(Y, length)

  ok <- ni > 0
  Y <- Y[ok]; bw.T <- as.list(qmid)[ok]

  Z <- mapply(density.default, x = Y, bw = bw.T, SIMPLIFY = F,
              MoreArgs = list(kernel = "gaussian", n = dimt,
                              from = range.t[1], to = range.t[2] ))
  zx <- Z[[1]]$x
  Edge0 <- sapply(bw.T, edge.t, dt = zx, TS = range.t[2], TI = range.t[1])
  Ys <- sapply(Z, "[[", 'y') / Edge0
  Ys <- apply(sweep(Ys, 2, ni[ok], FUN = "*"), 1, sum)
  fz.t <- approxfun(x = zx, y = Ys)

  ZZ <- switch(at,
               points = fz.t(t),
               bins = data.frame(x = zx, y = Ys))
  return(ZZ)
}

edge.t <- function(bwt, dt, TS = max(dt), TI = min(dt)) {
  pnorm((TS - dt) / bwt) - pnorm((TI- dt) / bwt)
}
