#'@importFrom spatstat.utils check.1.integer

#Check groups and things
check.ngroups <- function(Ngroups, N, order = 2){
  if (missing(Ngroups) || is.null(Ngroups)) {
    ng <- max(1L, floor(N ^ (1 / order)))
  }
  else if (any(is.infinite(Ngroups))) {
    ng <- N
  }
  else {
    check.1.integer(Ngroups)
    ng <- min(N, Ngroups)
  }
  return(ng)
}