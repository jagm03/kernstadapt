#' @importFrom sparr OS
#' @importFrom stats quantile
#' @importFrom spatstat.geom %mark% marks
#' @importFrom spatstat.random rpoispp
#' @export
dens.st.bin <- function(X, t,
                        bw.xy = NULL, bw.t = NULL, #bandwidths
                        ngroups.xy = NULL, ngroups.t = NULL, #groups
                        dimt = 128, dimyx = 128, #resolution
                        at = c("points", "bins") #at
){
  at <- match.arg(at)
  nT <- length(t)
  if(nT != X$n) stop(paste("Length of temporal vector does not match number of spatial observations"))
  range.t <- range(t)
  ngroups.xy <- check.ngroups(ngroups.xy, N = nT, order = 3)
  ngroups.t <- check.ngroups(ngroups.t, N = nT, order = 6)

  if (missing(bw.t) || is.null(bw.t)) {
    bw.t <- bw.abram.temp(t)
  }
  else if (is.numeric(bw.t)) {
    check.nvector(bw.t, nT, oneok = TRUE)
    if (length(bw.t) == 1)
      bw.t <- rep(bw.t, nT)
  }
  else stop("Argument 'bw.t' should be a single value or a numeric vector")

  if (missing(bw.xy) || is.null(bw.xy)) {
    bw.xy <- bw.abram(X, h0 = OS(X))
  }
  else if (is.numeric(bw.xy)) {
    check.nvector(bw.xy, nT, oneok = TRUE)
    if (length(bw.xy) == 1)
      bw.xy <- rep(bw.xy, nT)
  }
  else stop("Argument 'bw.xy' should be a single value or a numeric vector")

  p.xy <- seq(0, 1, length = ngroups.xy + 1)
  p.t <- seq(0, 1, length = ngroups.t + 1)
  qbands.xy <- quantile(bw.xy, p.xy)
  qbands.t <- quantile(bw.t, p.t)
  groupid.xy <- findInterval(bw.xy, qbands.xy, all.inside = T)
  groupid.t <- findInterval(bw.t, qbands.t, all.inside = T)
  pmid.xy <- (p.xy[ -1 ] + p.xy[ -length(p.xy) ]) / 2
  pmid.t <- (p.t[ -1 ] + p.t[ -length(p.t) ]) / 2
  qmid.xy <- quantile(bw.xy, pmid.xy)
  qmid.t <- quantile(bw.t, pmid.t)
  group.xy <- factor(groupid.xy, levels = 1:ngroups.xy)
  group.t <- factor(groupid.t, levels = 1:ngroups.t)
  BX <- X %mark% data.frame(t = t, g.t = group.t)
  Y.xy <- split(BX, group.xy)

  # remove zero pp
  ok1 <- sapply(Y.xy, npoints) > 0
  Y.xy <- Y.xy[ok1]
  qmid.xy <- qmid.xy[ok1]

  Y.xy <- mapply(function(x, e) {
    marks(x) <- cbind(marks(x), epsilon = e)
    return(x)
  },
  x = Y.xy, e = as.list(qmid.xy), SIMPLIFY = F)

  Y.xy.t <- unlist(lapply(Y.xy, split, f = "g.t"), recursive = F)
  Y.xy.t <- lapply(Y.xy.t, function(x) {
    x$marks$delta = qmid.t[x$marks$g.t]
    return(x)})

  PP <- lapply(Y.xy.t, function(x) unmark(x) %mark% x$marks$t)
  B.XY <- lapply(Y.xy.t, function(x) marks(x)$epsilon[1])
  B.t <- lapply(Y.xy.t, function(x) marks(x)$delta[1])

  ok <- sapply(PP, npoints) > 0
  PP <- PP[ok]; B.XY <- B.XY[ok]; B.t <- B.t[ok]

  Z <- mapply(dens.direct.engine, X = PP, bw.xy = B.XY, bw.t = B.t, SIMPLIFY = F,
              MoreArgs = list(tlim = range.t, sres = dimyx, tres = dimt))

  ZZ <- Reduce("+", Z)
  WM <- as.mask(X$window, dimyx = dimyx)
  inside <- WM$m
  grx <- WM$xcol
  gry <- WM$yrow
  zz <- list()

  for(i in 1:dimt){
    zz[[i]] <- im(ZZ[,, i], xcol = grx, yrow = gry)
    zz[[i]][!inside] <- NA
  }

  if (at == "bins")
    return(zz)

  if (at == "points"){
    lambda <- rep(NA, nT)
    tcw <- diff(range.t) / dimt
    grt <- range.t[1] + 0.5 * tcw + (0:(dimt - 1)) * tcw
    tC <- findInterval(t, grt, all.inside = T)
    for (i in 1:nT) {
      lambda[i] <- safelookup(zz[[tC[i]]], X[i])
    }
    return(lambda)
  }
}
