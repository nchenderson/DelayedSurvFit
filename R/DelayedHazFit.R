DelayedHazFit <- function(times, events, trt, gamma.fixed=NULL, theta.fixed=NULL, max.times=100,
                           inner.iter=50, final.iter=1000, verbose=TRUE) {
  
  options(nloptr.show.inequality.warning=FALSE)
  num.unique <- length(unique(times))
  if(num.unique > max.times) {
    new.times <- DiscretizeTimes(times, max.times=max.times)
  } else {
    new.times <- times
  }
  sfit <- survfit(Surv(new.times, events) ~ trt)
  
  utimes <- unique(sort(new.times))
  risk.data <- data.frame(strata = summary(sfit, times = utimes, extend = TRUE)$strata, 
                          time = summary(sfit, times = utimes, extend = TRUE)$time, 
                          n.risk = summary(sfit, times = utimes, extend = TRUE)$n.risk,
                          n.event = summary(sfit, times = utimes, extend = TRUE)$n.event,
                          surv = summary(sfit, times=utimes, extend=TRUE)$surv)
  
  n1 <- risk.data$n.risk[risk.data$strata=="trt=1"]
  d1 <- risk.data$n.event[risk.data$strata=="trt=1"]
  S1 <- risk.data$surv[risk.data$strata=="trt=1"]
  n0 <- risk.data$n.risk[risk.data$strata=="trt=0"]
  d0 <- risk.data$n.event[risk.data$strata=="trt=0"]
  S0 <- risk.data$surv[risk.data$strata=="trt=0"]
  
  #max.index <- min(max(n0 > 0), max(n1 > 0))
  num.zeros0 <- sum(n0 - d0 == 0) 
  num.zeros1 <- sum(n1 - d1 == 0)
  if(num.zeros0 > 0 | num.zeros1 > 0) {
    if(num.zeros0 > 0 & num.zeros1 == 0) {
      cut.index <- min(which(n0 - d0 == 0))
    } else if(num.zeros0 == 0 & num.zeros1 > 0) {
      cut.index <- min(which(n1 - d1 == 0))
    } else if(num.zeros0 == 0 & num.zeros1 == 0) {
      cut.index <- min(min(which(n0 - d0==0)), min(which(n1 - d1==0)))
    }
    remove.indices <- cut.index:length(utimes)
    n1 <- n1[-remove.indices]
    d1 <- d1[-remove.indices]
    S1 <- S1[-remove.indices]
    
    n0 <- n0[-remove.indices]
    d0 <- d0[-remove.indices]
    S0 <- S0[-remove.indices]
    utimes <- utimes[-remove.indices]
  }
  
  n.pars <- 2*length(utimes)
  par.target <- c(log(n0 - d0) - log(n0), log(n1 - d1) - log(n1))
  
  ## need to clean data further here.
  
  nsqp1 <- inner.iter
  nsqp2 <- final.iter
  
  ## (3) Find optimal theta (if theta is not specified)
  if(is.null(theta.fixed) & is.null(gamma.fixed)) {
    theta.grid <- c(0, utimes)
    ngrid <- length(theta.grid)
    ell1 <- ellneg1 <- rep(0, ngrid)
    for(k in 1:ngrid) {
      ell1[k] <- SurvFnKnownTheta_Haz(theta = theta.grid[k], gamma=1, d0=d0, d1=d1, 
                                  n0 = n0, n1 = n1, utimes=utimes, max.sqp.iter=nsqp1)$loglik.val
      ellneg1[k] <- SurvFnKnownTheta_Haz(theta = theta.grid[k], gamma=-1, d0=d0, d1=d1, 
                                     n0 = n0, n1 = n1, utimes=utimes, max.sqp.iter=nsqp1)$loglik.val
      
       if(verbose) {
         cat(k, " out of ", ngrid, " iterations \n")
       }
    }
    ## Need to polish the solution somehow.
    
    best.gamma <- ifelse(min(ellneg1) < min(ell1), -1, 1)
    if(best.gamma == -1) {
      best.theta <- theta.grid[which.min(ellneg1)]
    } else if(best.gamma == 1) {
      best.theta <- theta.grid[which.min(ell1)]
    }
    ## If best.theta equals maximum value, just switch to 0
    if(best.theta==max(theta.grid)) {
      best.theta <- 0
      best.gamma <- ifelse(best.gamma == 1, -1, 1)
    }
  } else if(is.null(theta.fixed) & !is.null(gamma.fixed)) {
    theta.grid <- c(0, utimes)
    ngrid <- length(theta.grid)
    ell <- rep(0, ngrid)
    for(k in 1:ngrid) {
      ell[k] <- SurvFnKnownTheta_Haz(theta = theta.grid[k], gamma=gamma.fixed, d0=d0, d1=d1, 
                                 n0 = n0, n1 = n1, utimes=utimes)$loglik.val
      if(verbose) {
        cat(k, " out of ", ngrid, " iterations \n")
      }
    }
    ## Need to polish the solution somehow.
    
    best.theta <- theta.grid[which.min(ell)]
    best.gamma <- gamma
    if(best.theta==max(theta.grid)) {
      best.theta <- 0
      best.gamma <- ifelse(best.gamma == 1, -1, 1)
    }
  } else if(!is.null(theta.fixed) & is.null(gamma.fixed)) {
    obj1 <- SurvFnKnownTheta_Haz(theta=theta.fixed, gamma=1, d0=d0, d1=d1, n0=n0, n1=n1, utimes=utimes)$loglik.val
    objneg1 <- SurvFnKnownTheta_Haz(theta=theta.fixed, gamma=-1, d0=d0, d1=d1, n0=n0, n1=n1, utimes=utimes)$loglik.val 
    
    best.gamma <- ifelse(obj1 < objneg1, 1, -1)
    best.theta <- theta.fixed
  } else if(!is.null(theta.fixed) & !is.null(gamma.fixed)) {
    theta.possible <- NULL
    best.theta <- theta.fixed
    best.gamma <- gamma.fixed
    objfn.check1 <- NULL
  }
  
  ## (4) 
  tmp <- SurvFnKnownTheta_Haz(theta=best.theta, gamma=best.gamma, d0=d0, 
                          d1=d1, n0=n0, n1=n1, utimes=utimes, max.sqp.iter=nsqp2) 
  
  ## Need to also return the discrete hazards themselves.
  ## note that this answer excludes the last jump point
  ans <- list(times=utimes, surv0=tmp$Surv0, surv1=tmp$Surv1, haz0=tmp$haz0, haz1=tmp$haz1,
              nevents0=d0, nrisk0=n0, nevents1=d1, nrisk1=n1, theta=best.theta, gamma=best.gamma, 
              discretized.times=new.times)
  class(ans) <- "surv.delay"
  return(ans)
}