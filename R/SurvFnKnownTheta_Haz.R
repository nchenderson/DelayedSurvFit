SurvFnKnownTheta_Haz <- function(theta, gamma, d0, d1, n0, n1, utimes, max.sqp.iter=100) {
  
  n.pars <- 2*length(utimes)
  
  Cmat <- ConstructConstrHazMat(utimes, theta=theta, gamma=gamma)
  Amat <- (-1)*rbind(Cmat, diag(1, n.pars))
  par.target <- c(log(n0 - d0) - log(n0), log(n1 - d1) - log(n1))
  
  qp.sol <- dykstra(Dmat=diag(rep(1, n.pars)), dvec=par.target, Amat=t(Amat))
  
  
  loglik.val.qp <- LogLik(qp.sol$solution, nevents0=d0, nevents1=d1, nrisk0=n0, nrisk1=n1) 
  
  Ashort <- (-1)*Cmat
  locconstrfn <- function(x) {
    #return(as.numeric(Amat%*%x))
    return(as.numeric(Ashort%*%x))
  }
  locconstrfnder <- function(x) {
    #return(Amat)
    return(Ashort)
  }
  sqp.sol <- slsqp(x0=qp.sol$solution, fn=LogLik, gr = LogLikDer, upper = rep(0, n.pars), hin = locconstrfn,
                   hinjac = locconstrfnder, control = list(maxeval=max.sqp.iter), nevents0=d0, nevents1=d1, nrisk0=n0, nrisk1=n1)
  loglik.val.sqp <- LogLik(sqp.sol$par, nevents0=d0, nevents1=d1, nrisk0=n0, nrisk1=n1) 
  
  if(loglik.val.qp <= loglik.val.sqp) {
    CS0 <- exp(cumsum(qp.sol$solution[1:(n.pars/2)]))
    CS1 <- exp(cumsum(qp.sol$solution[(n.pars/2 + 1):n.pars]))
    haz0 <- 1 - exp(qp.sol$solution[1:(n.pars/2)])
    haz1 <- 1 - exp(qp.sol$solution[(n.pars/2 + 1):n.pars])
  } else if(loglik.val.qp > loglik.val.sqp) {
   # print(summary(sqp.sol$par))
    
    CS0 <- exp(cumsum(sqp.sol$par[1:(n.pars/2)]))
    CS1 <- exp(cumsum(sqp.sol$par[(n.pars/2 + 1):n.pars]))
    haz0 <- 1 - exp(sqp.sol$par[1:(n.pars/2)])
    haz1 <- 1 - exp(sqp.sol$par[(n.pars/2 + 1):n.pars])
  }    
  ans <- list(Surv0 = CS0, Surv1 = CS1, haz0=haz0, haz1=haz1,
              loglik.val=loglik.val.qp)
  return(ans)
}
