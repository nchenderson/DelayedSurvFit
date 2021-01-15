rmst <- function(obj, tau0=NULL, tau1=NULL)
{
  if(is.null(tau0)) {
      tau0 <- max(obj$times)
  }
  if(is.null(tau1)) {
      tau1 <- max(obj$times)
  }
  ind0 <- obj$times <= tau0
  ind1 <- obj$times <= tau1
  t0 <- c(0, obj$times[ind0], tau0)
  t1 <- c(0, obj$times[ind1], tau1)
  
  S0 <- c(1, obj$surv0[ind0])
  S1 <- c(1, obj$surv1[ind1])
  
  rmst0 <- sum(S0*diff(t0))
  rmst1 <- sum(S1*diff(t1))
  
  return(list(rmst0=rmst0, rmst1=rmst1))
}


