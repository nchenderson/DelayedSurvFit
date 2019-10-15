ProfileLogLikSQP <- function(theta, nevents0, nevents1, nrisk0, nrisk1, utimes0, utimes1,
                             H0, H1) {
  utimes0.orig <- utimes0
  utimes1.orig <- utimes1
  
  a.set <- FindActiveSet(theta=theta, utimes0, utimes1) 
  nevents0 <- nevents0[a.set$active.set0]
  nrisk0 <- nrisk0[a.set$active.set0]
  utimes0 <- utimes0[a.set$active.set0]
  
  nevents1 <- nevents1[a.set$active.set1]
  nrisk1 <- nrisk1[a.set$active.set1]
  utimes1 <- utimes1[a.set$active.set1]
  
  
  na.w0 <- diff(H0(c(0, utimes0)))
  na.w1 <- diff(H1(c(0, utimes1)))
  d.target <- c(na.w0, na.w1)
  
  Cmat <- ConstructConstrMat(utimes0, utimes1, theta)
  n.pars <- length(utimes0) + length(utimes1)
  
  Dmat <- rbind(Cmat, diag(rep(1, n.pars)), diag(rep(-1, n.pars)))
  bvec <- c(rep(0, nrow(Cmat)), rep(0, n.pars), rep(-1, n.pars))
  Amat <- t(Dmat)
  a <- solve.QP(diag(rep(1, n.pars)), dvec=d.target, Amat=Amat, bvec=bvec)
  par.init <- a$solution
  par.init <- par.init + 1e-12 ## sometimes solve.QP returns very small negative numbers (e.g., -1e-19)
    ### All the w0 between theta and min(utimes1: utimes1 > theta) should also be zero.
  
  ## need to ensure that this initial value is feasible.
  niter <- 200
  old.par <- par.init
  Dmat <- diag(old.par)
  resid.sq <- 1
  k <- 1
  loglik.old <- LogEL(par.init, nevents0, nevents1, nrisk0, nrisk1) 
  print(loglik.old)
  alpha <- 0.5
  ## Use SQP for finding optimal solution.
  while(resid.sq > 1e-6 & k <= niter) {
    dd <- LogELDer2(old.par, nevents0, nevents1, nrisk0, nrisk1) 
    R <- LogELDer(old.par, nevents0, nevents1, nrisk0, nrisk1)
    dvec <- (-1)*R
    diag(Dmat) <- dd
    # Is this the right quadratic approximation here?
    
    a <- solve.QP(Dmat, dvec, Amat, bvec)
    
    #a <- MyQP(dd, qvec=-dvec, Amat=-t(Amat), bvec=bvec)
    
    new.par <- old.par + alpha*a$solution
    #new.par <- old.par + alpha*a
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
  
  cum.H0 <- stepfun(utimes0, cumsum(c(0,hazard0)), right=FALSE)
  cum.H1 <- stepfun(utimes1, cumsum(c(0,hazard1)), right=FALSE)
  
  tt0 <- length(utimes0.orig)
  tt1 <- length(utimes1.orig)
  ff0 <- (H0(utimes0.orig[-tt0]) - cum.H0(utimes0.orig[-tt0]))^2
  ff1 <- (H1(utimes1.orig[-tt1]) - cum.H1(utimes1.orig[-tt1]))^2
  
  I0 <- sum(diff(utimes0.orig)*ff0)
  I1 <- sum(diff(utimes1.orig)*ff1)
  
  DistNA <- I0 + I1
  #ans <- LogEL(new.par, nevents0, nevents1, nrisk0, nrisk1) 
  return(DistNA)
}