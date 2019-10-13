
LogEL <- function(par, nevents0, nevents1, nrisk0, nrisk1, ...) {
  #eps <- 0
  par <- par 
  n0 <- length(nevents0)
  n1 <- length(nevents1)
  w0 <- par[1:n0] 
  w1 <- par[(n0 + 1):(n0 + n1)]
  ans0 <- sum(nevents0*log(w0)) + sum((nrisk0 - nevents0)*log(1 - w0))
  ans1 <- sum(nevents1*log(w1)) + sum((nrisk1 - nevents1)*log(1 - w1))
  return((-1)*ans0 - ans1)
}

LogELDer <- function(par, nevents0, nevents1, nrisk0, nrisk1, ...) {
  #eps <- 0
  par <- par 
  ans0 <- c(nevents0, nevents1)/par - (c(nrisk0,nrisk1) - c(nevents0,nevents1))/(1 - par)
  return(-ans0)
}

LogELDer2 <- function(par, nevents0, nevents1, nrisk0, nrisk1) {
  #eps <- 0
  par <- par 
  ans0 <- c(nevents0, nevents1)/(par*par) + (c(nrisk0,nrisk1) - c(nevents0,nevents1))/((1 - par)^2)
  return(ans0)
}
