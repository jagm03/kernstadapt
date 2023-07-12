#' Abramson's adaptive bandwidth for networks
#'
#' Computes adaptive smoothing bandwidth in the network case according to the inverse-square-root rule of Abramson (1982).
#'
#' @param X A point pattern on a linear network (object of class "lpp").
#' @param h0 The global smoothing bandwidth. The default is the maximal oversmoothing principle of Terrell (1990).
#' @param trim A trimming value to cut extreme large bandwidths.
#' @param at Character string specifying whether to compute bandwidths at the points (at = "points", the default) or to compute bandwidths at every bin in a bin grid (at = "bins").
#' @param smoother Smoother for the pilot. A function or character string, specifying the function to be used to compute the pilot estimate when pilot is NULL or is a point pattern.
#' @param ... Additional arguments passed to smoother to control the type of smoothing.
#'
#' @details
#' This function returns a set of adaptive smoothing bandwidths driven by Abramson's (1982) method for a point pattern in a linear network.
#' The bandwidth at location \eqn{u} is given by
#' \deqn{
#' \epsilon(u) = h0 * \mbox{min}\left[ \frac{1}{\gamma} \sqrt{\frac{n}{\tilde{lambda}(u)}}, \mbox{\texttt{trim}} \right]
#' }{
#' h0 * pmin((pilotvalues ^ (-0.5)) / gamma, trim)
#' }
#' where \eqn{\tilde{\lambda}(u)} is a pilot estimate of the network varying intensity and \eqn{\gamma} is a scaling constant depending on the pilot estimate.
#'
#' @return
#' If \code{at = "points"} (the default), the result is a numeric vector with one bandwidth for each data point in \code{X}.
#' If \code{at = "pixels"}, the output is an object with class "linim"
#'
#' @references
#'Abramson, I. (1982) On bandwidth variation in kernel estimates --- a square root law.
#' \emph{Annals of Statistics}, \bold{10}(4), 1217-1223.\cr
#'
#' @author Jonatan A. González
#'
#' @examples
#' #To be done
#'
#' @importFrom spatstat.geom is.lpp compatible.im
#' @importFrom spatstat.explore resolve.2D.kernel
#' @importFrom spatstat.linnet as.linim flatdensityfunlpp integral.linim densityQuick.lpp eval.linim
#' @export
bw.abram.net <- function(X, h0,
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
  } else if(!is.null(pilot)){
    stop("if supplied, 'pilot' must be a pixel image on the linear network or a point pattern",
         call.=FALSE)}

  if(!is.linim(pilot)) {
    if(is.character(smoother)) {
      smoother <- get(smoother, mode = "function")
    } else stopifnot(is.function(smoother))
    pilot <- smoother(pilot.data, sigma = hp, positive = TRUE, ...)
  }

  pilot <- pilot / integral.linim(pilot) # scale to probability density
  pilotvalues <- safelookup(pilot, as.ppp(pilot.data), warn=FALSE)
  ## geometric mean re-scaler (Silverman, 1986; ch 5).
  gamma <- exp(mean(log(pilotvalues[pilotvalues > 0])))^(-0.5)

  switch(at,
         points = {
           pilot.X <- safelookup(pilot, as.ppp(X), warn = FALSE)
           bw <- h0 * pmin((pilot.X^(-0.5)) / gamma, trim)
         },
         pixels = {
           bw <- eval.linim(h0 * pmin((pilot^(-0.5)) / gamma, trim))
         })

  return(bw)
}
