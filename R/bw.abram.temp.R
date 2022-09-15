#' Abramson's adaptive temporal bandwidths
#'
#' Computes adaptive smoothing bandwidth in the temporal case according to the inverse-square-root rule of Abramson (1982).
#'
#' @param t A vector (a temporal point pattern) from which the bandwidths should be computed.
#' @param h0 The global smoothing bandwidth. The default is Silverman's rule of thumb (bw.nrd0).
#' @param nt The number of equally spaced points at which the temporal density is to be estimated.
#' @param trim A trimming value to cut extreme large bandwidths.
#' @param at Character string specifying whether to compute bandwidths at the points (at = "points", the default) or to compute bandwidths at every bin in a bin grid (at = "bins").
#'
#' @details
#' This function returns a set of temporal adaptive smoothing bandwidths driven by the methods of Abramson (1982) and Hall and Marron (1988).
#' The bandwidth at location \eqn{v} is given by
#' \deqn{
#' \delta(v) = h0 * \mbox{min}\left[ \frac{1}{\gamma} \sqrt{\frac{n}{\lambda^{\mbox{t}}(v)}}, \mbox{\texttt{trim}} \right]
#' }{
#' h0 * pmin((pilotvalues ^ (-0.5)) / gamma, trim)
#' }
#' where \eqn{\lambda^{\mbox{t}}(v)} is a pilot estimate of the temporally varying intensity and \eqn{\gamma} is a scaling constant depending on the pilot estimate.
#'
#' @return
#' If \code{at = "points"} (the default), the result is a numeric vector with one bandwidth for each data point in \code{X}.
#' If \code{at = "bins"}, the output is an object with class "density" where \code{y} component is a vector with the estimated intensity values (see \link[stats]{density}).
#'
#' @references
#'Abramson, I. (1982) On bandwidth variation in kernel estimates --- a square root law.
#' \emph{Annals of Statistics}, \bold{10}(4), 1217-1223.\cr
#'
#' Davies, T.M. and Baddeley, A. (2018) Fast computation of spatially adaptive kernel estimates.
#' \emph{Statistics and Computing}, \bold{28}(4), 937-956.\cr
#'
#' Davies, T.M., Marshall, J.C., and Hazelton, M.L. (2018)
#' Tutorial on kernel estimation of continuous spatial and spatiotemporal relative risk.
#' \emph{Statistics in Medicine}, \bold{37}(7), 1191-1221.\cr
#'
#' Hall, P. and Marron, J.S. (1988) Variable window width kernel density estimates of probability
#' densities. \emph{Probability Theory and Related Fields}, \bold{80}, 37-49.\cr
#'
#' Silverman, B.W. (1986) \emph{Density Estimation for Statistics and Data Analysis}.
#' Chapman and Hall, New York.
#'
#' González J.A. and Moraga P. (2018)
#' An adaptive kernel estimator for the intensity function of spatio-temporal point processes
#' <arXiv:2208.12026>

#' @author Jonatan A. González
#'
#' @examples
#' t <- 2 * rbeta(100, 1.5, 5.5, 0.2) #Simulated temporal point pattern
#' bw.abram.temp(t)
#'
#' @importFrom stats bw.nrd0 density.default approxfun
#' @importFrom spatstat.utils check.1.real
#' @export
bw.abram.temp <- function (t, h0 = NULL, nt = 128, trim = NULL, at = "points")
{
  t <- checkt(t)
  stopifnot(sum(t < 0) == 0)
  if (missing(trim) || is.null(trim)) {
    trim <- 0.25 * diff(range(t))
  }
  if (missing(h0) || is.null(h0)) {
    h0 <- bw.nrd0(t)
  }
  else {
    check.1.real(h0)
    stopifnot(h0 > 0)
  }
  check.1.real(trim)
  stopifnot(trim > 0)
  pilot.data <- t
  pilot <- density.default(t, bw = h0, kernel = "gaussian", n = nt)
  pilotfun <- approxfun(pilot$x, pilot$y)
  pilotvalues <- pilotfun(pilot.data)
  gamma <- exp(mean(log(pilotvalues[pilotvalues > 0]))) ^ (- 0.5)
  switch(at, points = {
    bw <- h0 * pmin((pilotvalues ^ (-0.5)) / gamma, trim)
  }, bins = {
    bw <- pilot
    bw$y <- (h0 * pmin((pilot$y ^ (-0.5)) / gamma, trim))
  })
  return(bw)
}
