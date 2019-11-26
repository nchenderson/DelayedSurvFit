logsumplus <- function(a, b) {
    ind <- a==0
    ans <- sum(a[!ind]*log(b[!ind]))
    return(ans)
}

LogEL <- function(par, nevents0, nevents1, nrisk0, nrisk1, ...) {
  #eps <- 0
  par <- par + 1e-12
  n0 <- length(nevents0)
  n1 <- length(nevents1)
  w0 <- par[1:n0] 
  w1 <- par[(n0 + 1):(n0 + n1)]
 # ans0 <- sum(nevents0*log(w0)) + sum((nrisk0 - nevents0)*log(1 - w0))
#  ans1 <- sum(nevents1*log(w1)) + sum((nrisk1 - nevents1)*log(1 - w1))
  ans0 <- sum(nevents0*log(w0)) + logsumplus(nrisk0 - nevents0, 1 - w0)
  ans1 <- sum(nevents1*log(w1)) + logsumplus(nrisk1 - nevents1, 1 - w1)
  return((-1)*ans0 - ans1)
}

LogELDer <- function(par, nevents0, nevents1, nrisk0, nrisk1, ...) {
  #eps <- 0
  par <- par + 1e-12
  RR <- c(nrisk0,nrisk1) - c(nevents0,nevents1)
  ind <- RR==0
  ans0 <- c(nevents0, nevents1)/par 
  ans1 <- RR
  ans1[ind] <- 0
  ans1[!ind] <- RR[!ind]/(1 - par[!ind])
  ans <- ans1 - ans0
  return(ans)
}

LogELDer2 <- function(par, nevents0, nevents1, nrisk0, nrisk1) {
  #eps <- 0
  par <- par + 1e-12
  RR <- c(nrisk0,nrisk1) - c(nevents0,nevents1)
  ind <- RR==0
  ans0 <- c(nevents0, nevents1)/(par*par) 
  ans1 <- RR
  ans1[ind] <- 0
  ans1[!ind] <- RR[!ind]/((1 - par[!ind])^2)
  ans <- ans0 + ans1
  return(ans)
}
