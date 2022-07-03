#' @importFrom spatstat.geom verifyclass cut.ppp quadratcount.ppp
X <- rpoispp(100)
t <- abs(rbeta(npoints(X), shape1 = 0.1, shape2 = 0.3))

separability.test(X, t, nx = NULL, ny = NULL, nt = NULL){
  verifyclass(X, "ppp")
  if (missing(nx) || is.null(nx)) {
    nx <- ceiling(npoints(X) ^ (1/6))
  }
  if (missing(nx) || is.null(nt)) {
    ny <- ceiling(npoints(X) ^ (1/6))
  }
  if (missing(nt) || is.null(nt)) {
    nt <- ceiling(npoints(X) ^ (1/6))
  }
  Tj <- split.ppp(cut(X, breaks = nt))
  nij <- sapply(Tj, function(a) as.vector(quadratcount.ppp(a, nx = nx, ny = ny)))


}
# Separability


fisher.test(nij, simulate.p.value = T,  B = 50000)
