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

checkt <- function(t){
  if(is.null(t)) stop("observation times must be supplied as 't' or as a numeric vector of marks to 'X'")
  if(!is.vector(t)) stop("'t' must be a vector")
  if(!is.numeric(t)) stop("'t' must be numeric")
  return(t)
}
