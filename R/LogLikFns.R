LogLik <- function(par, nevents0, nevents1, nrisk0, nrisk1, ...) {
  ## Note that this returns the value of the negative log-likelihood
  m <- length(nevents0)
  ind0 <- nevents0 != 0
  ind1 <- nevents1 != 0
  par0 <- par[1:m]
  par1 <- par[(m+1):(2*m)]
  eps <- 1e-12
  
  t0 <- sum(nevents0[ind0]*log(-expm1(par0[ind0]) + eps)) + sum((nrisk0 - nevents0)*par0)
  t1 <- sum(nevents1[ind1]*log(-expm1(par1[ind1]) + eps)) + sum((nrisk1 - nevents1)*par1)
  return(-t0 - t1)
}


LogLikDer <- function(par, nevents0, nevents1, nrisk0, nrisk1, ...) {
  #eps <- 0
  m <- length(nevents0)
  par0 <- par[1:m]
  par1 <- par[(m+1):(2*m)]
  t0 <- t1 <- rep(0, m)
  
  ind0 <- nevents0 !=0
  ind1 <- nevents1 != 0
  t0[ind0] <- -nevents0[ind0]*(exp(par0[ind0])/(1 - exp(par0[ind0]))) + (nrisk0[ind0] - nevents0[ind0])
  t0[!ind0] <-  (nrisk0[!ind0] - nevents0[!ind0])
  t1[ind1] <- -nevents1[ind1]*(exp(par1[ind1])/(1 - exp(par1[ind1]))) + (nrisk1[ind1] - nevents1[ind1])
  t1[!ind1] <-  (nrisk1[!ind1] - nevents1[!ind1])
  #t1 <- -nevents1*expit(par1) + (nrisk1 - nevents1)
  ans <- c(-t0, -t1)
  return(ans)
}



