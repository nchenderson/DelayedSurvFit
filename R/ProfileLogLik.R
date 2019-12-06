ProfileLogLikSQP <- function(theta, gamma, nevents0, nevents1, nrisk0, nrisk1, utimes0, utimes1,
                             H0, H1) {
  
  ################################################################
  ## (1) Compute active set and constraint matrix
  utimes0.orig <- utimes0
  utimes1.orig <- utimes1
  
  a.set <- FindActiveSet(theta=theta, gamma, utimes0, utimes1) 
  nevents0 <- nevents0[a.set$active.set0]
  nrisk0 <- nrisk0[a.set$active.set0]
  utimes0 <- utimes0[a.set$active.set0]
  
  nevents1 <- nevents1[a.set$active.set1]
  nrisk1 <- nrisk1[a.set$active.set1]
  utimes1 <- utimes1[a.set$active.set1]
  
  na.w0 <- diff(H0(c(0, utimes0)))
  na.w1 <- diff(H1(c(0, utimes1)))
  d.target <- c(na.w0, na.w1)
  
  Cmat <- ConstructConstrMat(utimes0, utimes1, theta, gamma)
  n.pars <- length(utimes0) + length(utimes1)
  
  ################################################################
  ## (2) Perform parameter initialization
  ##       The constraints are Cw >= 0, w >= 0, w < 1
  Dmat <- rbind(Cmat, diag(rep(1, n.pars)), diag(rep(-1, n.pars)))
  bvec <- c(rep(0, nrow(Cmat)), rep(0, n.pars), rep(-1, n.pars))
  Amat <- t(Dmat)
  
  #a <- solve.QP(diag(rep(1, n.pars)), dvec=d.target, Amat=Amat, bvec=bvec)
  PP <- sparseMatrix(i=1:n.pars, j=1:n.pars, x=rep(1, n.pars), dims=c(n.pars,n.pars))
  
  a <- solve_osqp(P = PP, q = -d.target, A = t(Amat), l = bvec, u = NULL, 
                  pars=osqpSettings(verbose=FALSE,eps_prim_inf = 1e-10))
  par.init1 <- a$x
  par.init1[par.init1 < 0] <- 1e-12
  a <- try(solve.QP(diag(rep(1, n.pars)), dvec=d.target, Amat=Amat, bvec=bvec))
  if(class(a)!='try-error') {
    par.init0 <- a$solution
    par.init0[par.init0 < 0] <- 1e-12
  } else {
    par.init0 <- par.init1
  }
  loglik0 <- LogEL(par.init0, nevents0, nevents1, nrisk0, nrisk1) 
  loglik1 <- LogEL(par.init1, nevents0, nevents1, nrisk0, nrisk1)
  if(!is.na(loglik1) & !is.na(loglik0)) {
     if(loglik1 < loglik0) {
        par.init <- par.init1
     } else {
        par.init <- par.init0
     }
  } else {
      par.init <- par.init0
  }
  #par.init <- pmax(par.init + 1e-12, 1e-12)
  par.init[par.init < 0] <- 1e-12
  
  tau0 <- length(nevents0)
  tau1 <- length(nevents1)
  niter <- 200
  old.par <- par.init
  Dmat <- diag(old.par)
  resid.sq <- 1
  k <- 1
  loglik.old <- LogEL(par.init, nevents0, nevents1, nrisk0, nrisk1) 
  alpha <- 0.5
  ## Use SQP for finding optimal solution.
  npars <- length(par.init)
  Dmat <- sparseMatrix(i=1:npars, j=1:npars, x=d.target, dims=c(npars,npars))
  while(resid.sq > 1e-6 & k <= niter) {
    dd <- LogELDer2(old.par, nevents0, nevents1, nrisk0, nrisk1) 
    R <- LogELDer(old.par, nevents0, nevents1, nrisk0, nrisk1)
    dvec <- (-1)*R
    diag(Dmat) <- dd
    # Is this the right quadratic approximation here?
    
    #a <- try(solve.QP(Dmat, dvec, Amat, bvec))
    a <- solve_osqp(P = Dmat, q = -dvec, A = t(Amat), l = bvec, u = NULL, pars=osqpSettings(verbose=FALSE))
    if(class(a)=="try-error") {
        a.tmp <- solve_osqp(P = Dmat, q = -dvec, A = t(Amat), l = bvec, u = NULL, pars=osqpSettings(verbose=FALSE))
        qp.solution <- a.tmp$x 
    } else {
        #qp.solution <- a$solution
        qp.solution <- a$x
    }

    new.par <- rep(-1, length(qp.solution))
    while(sum(new.par < -1e-12) > 0 | sum(new.par > 1) > 0) {
       new.par <- pmin(old.par + alpha*qp.solution, 1)
       alpha <- alpha/2
    }
    loglik.new <- LogEL(new.par, nevents0, nevents1, nrisk0, nrisk1)  
    if(loglik.new < loglik.old) {
      alpha <- 0.5
    } else {
      alpha <- alpha/2
      new.par <- pmin(old.par + alpha*qp.solution, 1)
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
  return(DistNA)
}