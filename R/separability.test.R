#' Separability test for spatio-temporal point processes
#'
#' Performs a separability test of the first-order intensity function based on a Fisher Monte Carlo test of cell counts.
#'
#' @param X A spatial point pattern (an object of class \code{ppp}) with the spatial coordinates of the observations.
#' @param t A numeric vector of temporal coordinates with equal length to the number of points in \code{X}. This gives the time associated with each spatial point.
#' @param nx,ny,nt Numbers of quadrats in the \eqn{x,y} and \eqn{t} directions.
#' @param nperm An integer specifying the number of replicates used in the Monte Carlo test.
#'
#' @details
#' This function performs a basic test of the separability hypothesis in a manner similar to independence test in two-way contingency tables.
#' The test is conditional on the observed number of points.
#' It considers a regular division of the interval \eqn{T} into disjoint sub-intervals \eqn{T_1,...,T_{n_t}} and similarly a division of the window \eqn{W} into disjoint subsets \eqn{W_1, ..., W{n_x \times n_y}}.
#' Then the function computes Fisher's test statistic and get a p-value based on Monte Carlos approximation.
#'
#' @return
#' A list with class "htest" containing the following components:
#' \item{p.value}{the approximate p-value of the test.}
#' \item{method}{the character string "Separability test based on Fisher's for counting data".}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{data.name}{a character string giving the name(s) of the data.}
#'
#' @references
#' Ghorbani et al. (2021) Testing the first-order separability hypothesis for spatio-temporal point patterns,
#' \emph{Computational Statistics & Data Analysis,} \bold{161}, p.107245.
#'
#' @author Jonatan A. González
#'
#' @note This is a fast preliminary separability test.
#'
#' @examples
#' data(lGCpp)
#' separability.test(lGCpp, nx = 5, ny = 4, nt = 3, nperm = 500)
#'
#' \donttest{
#' data(aegiss)
#' separability.test(aegiss, nx = 8, ny = 8, nt = 4)
#' }
#'
#' @importFrom spatstat.geom verifyclass cut.ppp quadratcount.ppp
#' @importFrom stats fisher.test
#' @export
separability.test <- function(X, t = NULL,
                              nx = NULL, ny = NULL, nt = NULL,
                              nperm = 1000)
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



