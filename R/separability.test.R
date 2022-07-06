#' @importFrom spatstat.geom verifyclass cut.ppp quadratcount.ppp
#' @importFrom stats fisher.test
#' @export
separability.test <- function(X, t = NULL, nx = NULL, ny = NULL, nt = NULL, nperm = 1000)
{
  verifyclass(X, "ppp")
  n <- npoints(X)
  if (missing(nx) || is.null(nx)) {
    nx <- ceiling(n ^ (1/6))
  }
  if (missing(nx) || is.null(nt)) {
    ny <- ceiling(n ^ (1/6))
  }
  if (missing(nt) || is.null(nt)) {
    nt <- ceiling(n ^ (1/6))
  }
  if(is.null(t)) t <- marks(X)
  t <- checkt(t)
  if(length(t) != n)
    stop(paste("Length of temporal vector does not match number of spatial observations\n   npoints(X) = ",n,"; length(t) = ",length(t), sep = ""))

  X <- unmark(X)
  X <- setmarks(X, t)
  Tj <- split.ppp(cut(X, breaks = nt))
  ok <- sapply(Tj, npoints) > 0   # removing zero pp
  Tj <- Tj[ok]
  nij <- sapply(Tj, function(a) as.vector(quadratcount.ppp(a, nx = nx, ny = ny)))
  testsep <- fisher.test(nij, simulate.p.value = T,  B = nperm)
  testsep$method[1] <- "Separability test based on Fisher's for counting data"
  testsep$alternative <- "Not spatio-temporal separability"
  testsep$data.name <- "Point pattern X with times t"
  testsep
}



