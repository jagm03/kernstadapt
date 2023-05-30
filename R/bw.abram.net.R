#' Abramson's adaptive bandwidth for networks
#'
#' Computes adaptive smoothing bandwidth in the network case according to the inverse-square-root rule of Abramson (1982).
#'
#' @param X A point pattern on a linear network (object of class "lpp").
#' @param h0 The global smoothing bandwidth. The default is the maximal oversmoothing principle of Terrell (1990).
#' @param trim A trimming value to cut extreme large bandwidths.
#' @param at Character string specifying whether to compute bandwidths at the points (at = "points", the default) or to compute bandwidths at every bin in a bin grid (at = "bins").
#' @param smoother Smoother for the pilot. A function or character string, specifying the function to be used to compute the pilot estimate when pilot is NULL or is a point pattern.
#' @param ... Aditional arguments passed to smoother to control the type of smoothing.
#'
#' @details
#' This function returns a set of adaptive smoothing bandwidths driven by Abramson's (1982) method.
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
#' González J.A. and Moraga P. (2018)
#' An adaptive kernel estimator for the intensity function of spatio-temporal point processes
#' <https://arxiv.org/pdf/2208.12026.pdf>

#' @author Jonatan A. González
#'
#' @examples
#' t <- 2 * rbeta(100, 1.5, 5.5, 0.2) #Simulated temporal point pattern
#' bw.abram.temp(t)
#'
#' @importFrom stats bw.nrd0 density.default approxfun
#' @importFrom spatstat.utils check.1.real
#' @importFrom spatstat.geom is.lpp
#' @importFrom spatstat.explore resolve.2D.kernel
#' @importFrom spatstat.linet integral.linim densityQuick.lpp
#' @export
bw.abram.ppp <- function(X, h0,
                         at = c("points", "pixels"),
                         hp = h0, pilot = NULL, trim = 5,
                         smoother = densityQuick.lpp,
                         ...){
  stopifnot(is.lpp(X))#
  at <- match.arg(at)

  if(missing(h0) || is.null(h0)) {
    h0 <- OS(as.ppp(X))
  } else {
    check.1.real(h0)
    stopifnot(h0 > 0)
  }

  check.1.real(trim)
  stopifnot(trim > 0)

  pilot.data <- X
  imwin <- as.linim(flatdensityfunlpp(X))

  if(is.linim(pilot)){
    if(!compatible.im(imwin, pilot))
      stop("'X' and 'pilot' have incompatible network domains", call.=FALSE)
    #' clip the worst small values away
    pilot[pilot <= 0] <- min(pilot[pilot>0])
  } else if(is.ppp(pilot)){
    if(!compatible.im(imwin, as.linim(flatdensityfunlpp(pilot))))
      stop("'X' and 'pilot' have incompatible network domains", call.=FALSE)
    pilot.data <- pilot
  } else if(!is.null(pilot))
    stop("if supplied, 'pilot' must be a pixel image on the linear network or a point pattern",
         call.=FALSE)

  if(!is.linim(pilot)) {
    if(is.character(smoother)) {
      smoother <- get(smoother, mode = "function")
    } else stopifnot(is.function(smoother))
    pilot <- smoother(pilot.data, sigma = hp, positive = TRUE, ...)
  }

  pilot <- pilot / integral.linim(pilot) # scale to probability density
  pilotvalues <- safelookup(pilot, pilot.data, warn=FALSE)
  ## geometric mean re-scaler (Silverman, 1986; ch 5).
  gamma <- exp(mean(log(pilotvalues[pilotvalues > 0])))^(-0.5)

  switch(at,
         points = {
           pilot.X <- safelookup(pilot,X,warn=FALSE)
           bw <- h0 * pmin((pilot.X^(-0.5))/gamma,trim)
         },
         pixels = {
           bw <- eval.im(h0 * pmin((pilot^(-0.5))/gamma, trim))
         })

  return(bw)
}
