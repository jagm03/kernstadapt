#' @importFrom spatstat.utils check.nvector
#' @importFrom spatstat.core bw.abram density.ppp
#' @importFrom spatstat.geom split.ppp im as.mask safelookup marks unmark Window npoints
#' @importFrom misc3d kde3d
#' @importFrom sparr OS
#' @importFrom stats pnorm
#' @export
dens.adapt.direct <- function(X, t,
                              bw.xy = NULL, bw.t = NULL, #bandwidths
                              dimt = 128, dimyx = 128, #resolution
                              at = c("points", "bins") #at
){
  at <- match.arg(at)
  nT <- length(t)
  if(nT != X$n) stop(paste("Length of temporal vector does not match number of spatial observations"))
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
  WM <- as.mask(X$window)
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
