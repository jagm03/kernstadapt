
dens.sep.direct<- function(X, t = NULL,
                           dimt = 128, dimyx = 128, #resolution
                           bw.xy = NULL, bw.t = NULL, #bandwidths
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
  range.t <- range(t)
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
  #Temporal part
  PPt <- split(t, f = factor(1:X$n))
  Z <- mapply(density.default, x = PPt, bw = bw.t, SIMPLIFY = F,
              MoreArgs = list(kernel = "gaussian", n = dimt,
                              from = range.t[1], to = range.t[2] ))
  zx <- Z[[1]]$x
  Edge0 <- sapply(bw.t, edge.t, dt = zx, TS = range.t[2], TI = range.t[1])
  Ys <- sapply(Z, "[[", 'y') / Edge0
  Ys <- apply(Ys, 1, sum)
  Mtemporal <- switch(at,
                      points = approxfun(x = zx, y = Ys)(t),
                      bins = Ys)

  #Spatial part
  at.s <- switch (at, bins = "pixels", points = "points")
  Mspatial <- densityAdaptiveKernel.ppp(unmark(X), bw = bw.xy, at = at.s,
                                        dimyx = dimyx, edge = T, ngroups = Inf)
  Mst <- switch(at,
                points = Mtemporal * Mspatial / length(t),
                bins = lapply(Mtemporal, function(x) eval.im(Mspatial * x / length(t)))
  )
  return(Mst)
}
