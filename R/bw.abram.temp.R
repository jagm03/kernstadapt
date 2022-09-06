#' Computes adaptive smoothing bandwidth in the temporal case according to the inverse-square-root rule of Abramson (1982).
#'
#'@param ti A vector (a temporal point pattern) for which the bandwidths should be computed.
#'@param h0 The global smoothing bandwidth. The default is Silverman's rule of thumb (bw.nrd0).
#'@param trim A trimming value to cut extreme large bandwidths.
#'@param nt The number of equally spaced points at which the temporal density is to be estimated.
#'@param at Character string specifying whether to compute bandwidths at the points (at = "points", the default) or to compute bandwidths at every bin in a bin grid (at = "bins").
#'@examples
#'X <- rbeta(100, 1.5, 5.5, 0.2) * 2 #Simulated temporal point pattern
#'bw.abram.temp(X)
#'@importFrom stats bw.nrd0 density.default approxfun
#'@importFrom spatstat.utils check.1.real
#'@export
bw.abram.temp <- function (ti, h0 = NULL, trim = 5, nt = 128, at = "points")
{
  ti <- checkt(ti)
  stopifnot(sum(ti < 0) == 0)
  if (missing(h0) || is.null(h0)) {
    h0 <- bw.nrd0(ti)
  }
  else {
    check.1.real(h0)
    stopifnot(h0 > 0)
  }
  check.1.real(trim)
  stopifnot(trim > 0)
  pilot.data <- ti
  pilot <- density.default(ti, bw = h0, kernel = "gaussian", n = nt)
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
