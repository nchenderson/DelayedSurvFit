CumHazKnownTheta <- function(theta, gamma, nevents0, nevents1, nrisk0, nrisk1, utimes0, utimes1) {
  
  a.set <- FindActiveSet(theta=theta, gamma=gamma, utimes0, utimes1) 
  nevents0 <- nevents0[a.set$active.set0]
  nrisk0 <- nrisk0[a.set$active.set0]
  utimes0 <- utimes0[a.set$active.set0]
  
  nevents1 <- nevents1[a.set$active.set1]
  nrisk1 <- nrisk1[a.set$active.set1]
  utimes1 <- utimes1[a.set$active.set1]
  
  na.w0 <- diff(H0(c(0, utimes0)))
  na.w1 <- diff(H1(c(0, utimes1)))
  d.target <- c(na.w0, na.w1)
  
  init.find <- FindInitialVectors(theta=theta, utimes0, utimes1) 
  
  #par.init <- c(init.find$w0, init.find$w1)
  #mm <- max(par.init)
  #while(mm > 1) {
  #  par.init <- par.init/2
  #  mm <- max(par.init)
  #}
  
  Cmat <- ConstructConstrMat(utimes0, utimes1, theta, gamma)
  n.pars <- length(utimes0) + length(utimes1)
  
  Dmat <- rbind(Cmat, diag(rep(1, n.pars)), diag(rep(-1, n.pars)))
  bvec <- c(rep(0, nrow(Cmat)), rep(0, n.pars), rep(-1, n.pars))
  Amat <- t(Dmat)
  a <- solve.QP(diag(rep(1, n.pars)), dvec=d.target, Amat=Amat, bvec=bvec)
  par.init <- a$solution + 1e-12
  
  ### All the w0 between theta and min(utimes1: utimes1 > theta) should also be zero.
  niter <- 200
  old.par <- par.init
  Dmat <- diag(old.par)
  resid.sq <- 1
  k <- 1
  loglik.old <- LogEL(par.init, nevents0, nevents1, nrisk0, nrisk1)  
  alpha <- 0.5
  ## Use SQP for finding optimal solution.
  while(resid.sq > 1e-6 & k <= niter) {
    dd <- LogELDer2(old.par, nevents0, nevents1, nrisk0, nrisk1) 
    R <- LogELDer(old.par, nevents0, nevents1, nrisk0, nrisk1)
    dvec <- (-1)*R
    diag(Dmat) <- dd
    # Is this the right quadratic approximation here?
    
    a <- solve.QP(Dmat, dvec, Amat, bvec)
    
    new.par <- old.par + alpha*a$solution
    loglik.new <- LogEL(new.par, nevents0, nevents1, nrisk0, nrisk1)  
    print(c(loglik.new, loglik.old))
    if(loglik.new < loglik.old) {
      alpha <- 0.5
    } else {
      alpha <- alpha/2
      new.par <- old.par + alpha*a$solution
      print('reject')
      ## add better safeguards here later!
    }
    resid.sq <- sum((new.par - old.par)*(new.par - old.par))
    old.par <- new.par
    loglik.old <- loglik.new
    k <- k+1
  }
  tau0 <- length(nevents0)
  tau1 <- length(nevents1)
  hazard0 <- new.par[1:tau0]
  hazard1 <- new.par[(tau0 + 1):(tau0 + tau1)]
  return(list(hazard0=hazard0, hazard1=hazard1, nevents0=nevents0, nevents1=nevents1,
              utimes0=utimes0, utimes1=utimes1, nrisk0=nrisk0, nrisk1=nrisk1))
}
