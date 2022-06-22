#' @export
dens.t.bin <- function(ti, bw = NULL, ngroups = NULL, dimt = 128, at = c("points", "bins")){
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
