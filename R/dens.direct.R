#' Direct estimation of non-separable adaptive spatio-temporal intensity estimator
#'
#' Provides an adaptive-bandwidth kernel estimate for spatio-temporal point patterns in a non-separable fashion by calculating the classical estimator, i.e., the slowest estimation.
#'
#' @param X A spatial point pattern (an object of class \code{ppp}) with the spatial coordinates of the observations
#' @param t A numeric vector of temporal coordinates with equal length to the number of points in \code{X}. This gives the time associated with each spatial point.
#' @param dimt Temporal bin vector dimension. The default is 128.
#' @param dimyx Spatial pixel resolution. The default is 128 for each axes.
#' @param bw.xy Numeric vector of spatial smoothing bandwidths for each point in \code{X}. By default this is computed usign \link[spatstat.core]{bw.abram}, with \code{h0} given by \link[sparr]{OS}.
#' @param bw.t Numeric vector of temporal smoothing bandwidths for each point in \code{t}. By default this is computed usign \link{bw.abram.temp}.
#' @param at String specifying whether to estimate the intensity at a mesh (\code{at = "bins"}) or only at the points of \code{X} (\code{at = "points"}).
#'
#' @details
#' This function computes a non-separable spatio-temporal adaptive kernel estimate of the intensity. It starts from a planar point pattern \code{X} and a vector of times \code{t} and apply a non-separable kernel estimator for each of the points of \code{X}. The arguments \code{bw.xy} and \code{bw.t} specify the smoothing bandwidth vectors to be applied to each of the points in \code{X} and \code{t}. They should be a numeric vectors of bandwidths.
#'
#' @return
#' If \code{at = "points"} (the default), the result is a numeric vector with one entry for each data point in \code{X}. if \code{at = "bins"} is a list named (by time-point) list of pixel images (\link[spatstat.geom]{im} objects) corresponding to the joint spatio-temporal intensity over space at each discretised time bin.
#' @author Jonatan A. González
#'
#' @examples
#' \dontrun{
#' X <- rpoispp(400)
#' t <- rbeta(X$n, 1,4,0.8)
#' stIntensity <- dens.direct(X, t, dimt = 32, at = "bins")
#' plot(as.imlist(stIntensity[13:16]), main = 'Non-separable direct example')
#' }
#'
#' @importFrom spatstat.utils check.nvector
#' @importFrom spatstat.core bw.abram density.ppp
#' @importFrom spatstat.geom split.ppp im as.mask safelookup marks unmark Window npoints
#' @importFrom misc3d kde3d
#' @importFrom sparr OS
#' @importFrom stats pnorm
#' @export
dens.direct <- function(X, t = NULL,#point patterns
                        dimt = 128, dimyx = 128, #resolution
                        bw.xy = NULL, bw.t = NULL, #bandwidths
                        at = c("points", "bins") #at
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

  PP <- split.ppp(X, f = factor(1:X$n))

  Z <- mapply(dens.direct.engine, X = PP, t = as.list(t), bw.xy = bw.xy,
              bw.t = bw.t, SIMPLIFY = F, MoreArgs = list(tlim = range.t,
                                                         sres = dimyx,
                                                         tres = dimt))
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

dens.direct.engine <- function(X, t = NULL, bw.xy = NULL, bw.t = NULL,
                               tlim = NULL, sres = 128, tres = 64){
  if (is.null(t)){
    t <- marks(X)
    X <- unmark(X)
  }
  W <- Window(X)
  n <- npoints(X)

  WM <- as.mask(W, dimyx = sres)
  inside <- WM$m
  grx <- WM$xcol
  gry <- WM$yrow

  if(is.null(tlim)) tlim <- range(t)
  tcw <- diff(tlim) / tres
  grt <- tlim[1] + 0.5 * tcw + (0:(tres-1)) * tcw
  kt <- c(tlim[1] + 0.5 * tcw, tlim[2] - 0.5 * tcw)

  fhat <- kde3d(x = X$x, y = X$y, z = t, h = c(bw.xy, bw.xy, bw.t),
                n = c(sres, sres, tres), lims = c(range(grx), range(gry), kt))

  sz <- density.ppp(X, sigma = bw.xy, edge = T, dimyx = sres, spill = 1)
  sq <- im(matrix(1, sres, sres), xcol = grx, yrow = gry)
  sq <- sz$edg
  sq[sq > 1] <- 1
  sq[!inside] <- NA

  tq <- rep(1, tres)
  nearedge <- 1:tres
  wellinside <- which(grt > (tlim[1] + 4 * bw.t) & grt < (tlim[2] - 4 * bw.t))
  if(length(wellinside) > 0) nearedge <- nearedge[-wellinside]
  for(i in nearedge)
    tq[i] <- pnorm(tlim[2], mean = grt[i], sd = bw.t) - pnorm(tlim[1], mean = grt[i], sd = bw.t)

  z <- array(0, dim = dim(fhat$d))
  for(i in 1:tres){
    z[,,i] <- t(fhat$d[,,i]) / (sq$v * tq[i])
  }
  return(n * z)
}
