ProfileLogLikCO <- function(theta, nevents0, nevents1, nrisk0, nrisk1, utimes0, utimes1) {
  
  a.set <- FindActiveSet(theta=theta, utimes0, utimes1) 
  nevents0 <- nevents0[a.set$active.set0]
  nrisk0 <- nrisk0[a.set$active.set0]
  utimes0 <- utimes0[a.set$active.set0]
  
  nevents1 <- nevents1[a.set$active.set1]
  nrisk1 <- nrisk1[a.set$active.set1]
  utimes1 <- utimes1[a.set$active.set1]
  
  init.find <- FindInitialVectors(theta=theta, utimes0, utimes1) 
  
  
  par.init <- c(init.find$w0, init.find$w1)
  mm <- max(par.init)
  while(mm >= 1) {
    par.init <- par.init/2
    mm <- max(par.init)
  }
  
  Cmat <- ConstructConstrMat(utimes0, utimes1, theta)
  n.pars <- length(utimes0) + length(utimes1)
  
  Dmat <- rbind(Cmat, diag(rep(1, n.pars)), diag(rep(-1, n.pars)))
  bvec <- c(rep(0, nrow(Cmat)), rep(0, n.pars), rep(-1, n.pars))
  tmp2 <- as.vector(Dmat%*%par.init - bvec)
  print(summary(tmp2))
  
  ### All the w0 between theta and min(utimes1: utimes1 > theta) should also be zero.
  
  ## need to ensure that this initial value is feasible.
  tmp <- constrOptim(theta=par.init, f=LogEL, grad=LogELDer, ui=Dmat, ci=bvec,
                     nevents0=nevents0, nevents1=nevents1, nrisk0=nrisk0, nrisk1=nrisk1,
                     outer.eps=1e-7)
  new.par <- tmp$par

  #ans <- LogEL(new.par, nevents0, nevents1, nrisk0, nrisk1) 
  ans <- LogEL(new.par, nevents0, nevents1, nrisk0, nrisk1) 
  print(ans)
  return(ans)
}