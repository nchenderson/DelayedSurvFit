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
  
  
  #areas0 <- time.diff0*c(1, trun.surv0)
  
  #trunc.time0 <- sort(c(obj$times0[ind0], tau0))
  #trun.surv0 <- obj$surv0[ind0]
  #time.diff0 <- diff(c(0, trunc.time0))
  #areas0 <- time.diff0*c(1, trun.surv0)
  
  #ind1 <- obj$times1 <= tau1
  #trunc.time1 <- sort(c(obj$times1[ind1], tau1))
  #trun.surv1 <- obj$surv1[ind1]
  #time.diff1 <- diff(c(0, trunc.time1))
  #areas1 <- time.diff1*c(1, trun.surv1)
  
  #rmst0 <- sum(areas0)
  #rmst1 <- sum(areas1)
  return(list(rmst0=rmst0, rmst1=rmst1))
}


