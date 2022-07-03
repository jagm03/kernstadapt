#' @importFrom spatstat.geom verifyclass cut.ppp quadratcount.ppp
#' @export
separability.test <- function(X, t, nx = NULL, ny = NULL, nt = NULL, nperm = 1000)
{
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
  X <- setmarks(X, t)
  Tj <- split.ppp(cut(X, breaks = nt))
  nij <- sapply(Tj, function(a) as.vector(quadratcount.ppp(a, nx = nx, ny = ny)))
  testsep <- fisher.test(nij, simulate.p.value = T,  B = nperm)
  testsep$method[1] <- "Separability test based on Fisher's for counting data"
  testsep$alternative <- "Not spatio-temporal separability"
  testsep$data.name <- paste("Point pattern", (substitute(X)),
                             "with times",  (substitute(t)))
  testsep
}



