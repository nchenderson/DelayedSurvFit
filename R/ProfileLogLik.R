ProfileLogLikSQP <- function(theta, nevents0, nevents1, nrisk0, nrisk1, utimes0, utimes1) {
  
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
  while(mm > 1) {
      par.init <- par.init/2
      mm <- max(par.init)
  }
  
  Cmat <- ConstructConstrMat(utimes0, utimes1, theta)
  n.pars <- length(utimes0) + length(utimes1)
  
  Dmat <- rbind(Cmat, diag(rep(1, n.pars)), diag(rep(-1, n.pars)))
  bvec <- c(rep(0, nrow(Cmat)), rep(0, n.pars), rep(-1, n.pars))
 
  
  ### All the w0 between theta and min(utimes1: utimes1 > theta) should also be zero.
  
  Amat <- t(Dmat)
  ## need to ensure that this initial value is feasible.
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
    #Amat <- t(rbind(Cmat, diag(rep(1, length(par.init)))))
    #bvec <- as.numeric((-1)*crossprod(Amat, old.par))
    
    a <- solve.QP(Dmat, dvec, Amat, bvec)
    
    #a <- MyQP(dd, qvec=-dvec, Amat=-t(Amat), bvec=bvec)
    
    new.par <- old.par + alpha*a$solution
    #new.par <- old.par + alpha*a
    loglik.new <- LogEL(new.par, nevents0, nevents1, nrisk0, nrisk1)  
    print(summary(a$solution))
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
  #print(IsFeasible(new.par, Cmat))
  print(k)
  ans <- LogEL(new.par, nevents0, nevents1, nrisk0, nrisk1) 
  return(ans)
}