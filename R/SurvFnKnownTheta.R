SurvFnKnownTheta <- function(theta, gamma, d0, d1, n0, n1, utimes) {
  
    n.pars <- 2*length(utimes)
  
    Cmat <- ConstructConstrMat(utimes, theta=theta, gamma=gamma)
    Amat <- (-1)*rbind(Cmat, diag(1, n.pars))
    par.target <- c(log(n0 - d0) - log(n0), log(n1 - d1) - log(n1))

    qp.sol <- dykstra(Dmat=diag(rep(1, n.pars)), dvec=par.target, Amat=t(Amat))
    
   
    #eps <- 1e-4
    #qp.sol <- dykstra(Dmat=diag(rep(1, n.pars)), dvec=par.target, Amat=t(Amat), 
    #                  bvec = rep(eps, nrow(Amat)))
    ## might need to "polish" the solution here
    loglik.val <- LogLik(qp.sol$solution, nevents0=d0, nevents1=d1, nrisk0=n0, nrisk1=n1) 
      
    CS0 <- exp(cumsum(qp.sol$solution[1:(n.pars/2)]))
    CS1 <- exp(cumsum(qp.sol$solution[(n.pars/2 + 1):n.pars]))
    ans <- list(Surv0 = CS0, Surv1 = CS1, loglik.val=loglik.val)
    return(ans)
}
    
    
  