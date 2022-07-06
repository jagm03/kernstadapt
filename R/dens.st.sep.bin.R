#' @importFrom spatstat.core densityAdaptiveKernel.ppp
#' @importFrom spatstat.geom eval.im
#' @export
dens.st.sep.bin <- function(X, t = NULL,
                            bw.xy = NULL, bw.t = NULL, #bandwidths
                            ngroups.xy = NULL, ngroups.t = NULL, #groups
                            dimt = 128, dimyx = 128, #resolution
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
  if (is.null(bw.xy)) {
    h0 <- OS(X)
    bw.xy <- bw.abram(X, h0)
  }
  Mtemporal <- dens.t.bin(ti = t, bw = bw.t, ngroups = ngroups.t,
                          dimt = dimt, at = at)
  Mspatial <- densityAdaptiveKernel.ppp(unmark(X), bw = bw.xy, at = at.s,
                                        dimyx = dimyx, edge = T, ngroups = ngroups.xy)
  Mst <- switch(at,
                points = Mtemporal * Mspatial / length(t),
                bins = lapply(Mtemporal$y, function(x) eval.im(Mspatial * x / length(t)))
  )
  return(Mst)
}
