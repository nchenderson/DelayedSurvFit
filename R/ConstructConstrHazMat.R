ConstructConstrHazMat <- function(utimes, theta, gamma) {
  ### Function for constructing the main constraint matrix C
  #   we should have C * par \leq 0
  #   This is the constraint matrix under the assumption that
  #   the control arm initially "dominates" the active treatment arms 
  #ConstructConstrMat(utimes, theta, gamma)
  ## Assume that 0 <= \theta <= max event time
  times.tot <- utimes
  
  tau <- length(times.tot)
  tmp <- tmp2 <- matrix(0, nrow=tau, ncol=tau)
  tmp[lower.tri(tmp, diag=TRUE)] <- 1
  tmp2[lower.tri(tmp2, diag=TRUE)] <- 1
  
  if(min(times.tot) > theta) {
    ## all times are greater than theta
    Kstar <- 0
  } else {
    ## 
    Kstar <- max(which(times.tot <= theta))
  }
  
  
  if(gamma==1) {
    d1 <- c(rep(1, Kstar), rep(-1, tau - Kstar))
    d0 <- c(rep(-1, Kstar), rep(1, tau - Kstar))
    A <- diag(d1)
    B <- diag(d0)
  } else if(gamma == -1) {
    d1 <- c(rep(-1, Kstar), rep(1, tau - Kstar))
    d0 <- c(rep(1, Kstar), rep(-1, tau - Kstar))
    A <- diag(d1)
    B <- diag(d0)
  }
  ans <- cbind(A, B)
  return(ans)
}