DelayedSurvFit <- function(times, events, trt, gamma, theta.fixed=NULL) {
  a0 <- survfit(Surv(times[trt==0], events[trt==0]) ~ 1)
  a1 <- survfit(Surv(times[trt==1], events[trt==1]) ~ 1)
  
  tau0 <- length(a0$time)
  tau1 <- length(a1$time)
  utimes0 <- a0$time[-tau0]
  nevents0 <- a0$n.event[-tau0]
  nrisk0 <- a0$n.risk[-tau0]
  
  utimes1 <- a1$time[-tau1]
  nevents1 <- a1$n.event[-tau1]
  nrisk1 <- a1$n.risk[-tau1]
  
  if(is.null(theta.fixed)) {
    theta.interval <- c(min(c(utimes0, utimes1)), .95*max(c(utimes0, utimes1)) )
    ## probably need to give more thought to the best choice for theta.interval
    opt.gam1 <- optimize(f=ProfileLogLik, interval=theta.interval, nevents0=nevents0, nevents1=nevents1, 
                         nrisk0=nrisk0, nrisk1=nrisk1, utimes0=utimes0, utimes1=utimes1)
    best.theta <- opt.gam1$minimum
  } else {
    best.theta <- theta.fixed
  }
  tmp <- CumHazKnownTheta(best.theta, nevents0, nevents1, nrisk0, nrisk1, utimes0, utimes1)
  hazard0 <- tmp$hazard0
  hazard1 <- tmp$hazard1
  surv0 <- cumprod(1 - hazard0)
  surv1 <- cumprod(1 - hazard1)
  ## note that this answer excludes the last jump point
  ans <- list(times0=utimes0, nevents0=nevents0, nrisk0=nrisk0, times1=utimes1, nevents1=nevents1,
              nrisk1=nrisk1, hazard0=hazard0, hazard1=hazard1, surv0=surv0, surv1=surv1, theta=best.theta)
  class(ans) <- "surv.delay"
  return(ans)
}