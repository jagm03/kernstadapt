#' Adaptive kernel estimate of intensity of a temporal point pattern
#'
#' Estimates the intensity of a point process with only temporal dimension by applying an adaptive (variable bandwidth) Gaussian edge-corrected kernel smoothing.
#'
#' @param ti Temporal point pattern, a vector with observations.
#' @param bw Numeric vector of smoothing bandwidths for each point in ti. The default is to compute bandwidths using \link{bw.abram.temp}.
#' @param ngroups Number of groups into which the bandwidths should be partitioned and discretised. The default is the square root (rounded) of the number of points of \code{ti}.
#' @param dimt Bin vector dimension. The default is 128.
#' @param at String specifying whether to estimate the intensity at bins points (\code{at = "bins"}) or only at the points of ti (\code{at = "points"}).
#'
#' @details
#' This function computes a temporally-adaptive kernel estimate of the intensity from a one-dimensional point pattern ti using the partitioning technique of Davies and Baddeley (2018).
#' The argument bw specifies the smoothing bandwidths to be applied to each of the points in X. It should be a numeric vector of bandwidths.
#' Let the points of \eqn{ti} be \eqn{t_1, ..., t_n} and the corresponding bandwidths \eqn{\sigma_1,...,\sigma_n}, then the adaptive kernel estimate of intensity at a location \eqn{t} is
#' \deqn{\lambda(t) = \sum_{i=1}^n \frac{K(t,t_i; \sigma_i)}{c(t; \sigma_i)}}
#' where \eqn{K()} is the Gaussian smoothing kernel.
#' The method partition the range of bandwidths into ngroups intervals, correspondingly subdividing the points of the pattern \code{ti} into ngroups sub-patterns according to bandwidth, and applying fixed-bandwidth smoothing to each sub-pattern. Specifying \code{ngroups = 1} is the same as fixed-bandwidth smoothing with bandwidth \code{sigma = median(bw)}.
#'
#' @return
#' If \code{at = "points"} (the default), the result is a numeric vector with one entry for each data point in \code{ti}. If \code{at = "bins"} the result is a data.frame containing the \eqn{x,y} coordinates of the intensity function.
#'
#' @references
#' Davies, T.M. and Baddeley, A. (2018) Fast computation of spatially adaptive kernel estimates. Statistics and Computing, 28(4), 937-956.
#'
#' @author
#' Jonatan A. González
#'
#' @examples
#' \dontrun{
#' ti <- rbeta(100, 1,4,0.8)
#' dens.par.temp(ti)
#' }
#'
#' @export
dens.par.temp <- function(ti, bw = NULL, dimt = 128,
                       ngroups = NULL,
                       at = c("points", "bins")){
  stopifnot(sum(ti < 0) == 0)
  at <- match.arg(at)
  nT <- length(ti)
  range.t <- range(ti)
  if (missing(ngroups) || is.null(ngroups)) {
    ngroups <- max(1L, floor(sqrt(nT)))
  }
  else if (any(is.infinite(ngroups))) {
    ngroups <- nT
  }
  else {
    check.1.integer(ngroups)
    ngroups <- min(nT, ngroups)
  }

  if (missing(bw) || is.null(bw)) {
    bw <- bw.abram.temp(ti)
  }
  else if (is.numeric(bw)) {
    check.nvector(bw, nT, oneok = TRUE)
    if (length(bw) == 1)
      bw <- rep(bw, nT)
  }
  else stop("Argument 'bw' should be a single value or a numeric vector")

  p <- seq(0, 1, length = ngroups + 1)
  qbands <- quantile(bw, p)
  groupid <- findInterval(bw, qbands, all.inside = T)
  pmid <- (p[-1] + p[ -length(p) ]) / 2
  qmid <- quantile(bw, pmid)
  group <- factor(groupid, levels = 1:ngroups)
  Y <- split(ti, group, drop = F)
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
               points = fz.t(ti),
               bins = data.frame(x = zx, y = Ys))
  return(ZZ)
}

edge.t <- function(bwt, dt, TS = max(dt), TI = min(dt)) {
  pnorm((TS - dt) / bwt) - pnorm((TI- dt) / bwt)
}
