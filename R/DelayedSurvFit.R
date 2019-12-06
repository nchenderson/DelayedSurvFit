DelayedSurvFit <- function(times, events, trt, gamma=NULL, theta.fixed=NULL) {
  
  ### (1) Extract num.risk and num.events for each treatment arm
  a0 <- survfit(Surv(times[trt==0], events[trt==0]) ~ 1)
  a1 <- survfit(Surv(times[trt==1], events[trt==1]) ~ 1)
  
  ## look at only event times, 
  ## However, numerically there may be some benefit in keeping 
  ## the censoring times (think about it) 
  tau0 <- length(a0$time)
  tau1 <- length(a1$time)
  idx.event0 <- a0$n.event > 0
  idx.event1 <- a1$n.event > 0
  utimes0 <- a0$time[idx.event0]
  nevents0 <- a0$n.event[idx.event0]
  nrisk0 <- a0$n.risk[idx.event0]
  
  utimes1 <- a1$time[idx.event1]
  nevents1 <- a1$n.event[idx.event1]
  nrisk1 <- a1$n.risk[idx.event1]
  
  ## (2) Construct Nelson-Aalen estimates for each treatment arm
  h0 <- nevents0/nrisk0
  h1 <- nevents1/nrisk1
  H0 <- stepfun(utimes0, cumsum(c(0,h0)), right=FALSE)
  H1 <- stepfun(utimes1, cumsum(c(0,h1)), right=FALSE)
  h0.max <- max(utimes0)
  h1.max <- max(utimes1)
  
  ## (3) Find optimal theta (if theta is not specified)
  if(is.null(theta.fixed) & is.null(gamma)) {
    theta.interval <- c(min(c(utimes0, utimes1)), .95*max(c(utimes0, utimes1)) )
    ## probably need to give more thought to the best choice for theta.interval
    #opt.gam1 <- optimize(f=ProfileLogLikSQP, interval=theta.interval, gamma=1, nevents0=nevents0, nevents1=nevents1, 
    #                     nrisk0=nrisk0, nrisk1=nrisk1, utimes0=utimes0, utimes1=utimes1,
    #                     H0=H0, H1=H1)
    #opt.gamneg1 <- optimize(f=ProfileLogLikSQP, interval=theta.interval, gamma=-1, nevents0=nevents0, nevents1=nevents1, 
    #                     nrisk0=nrisk0, nrisk1=nrisk1, utimes0=utimes0, utimes1=utimes1,
    #                     H0=H0, H1=H1)
    theta.min <- min(c(utimes0, utimes1))/2
    theta.possible <- sort(c(theta.min, utimes0[-length(utimes0)], utimes1[-length(utimes1)]))

    if(length(theta.possible) > 200) {
        theta.possible <- seq(min(theta.possible), max(theta.possible), length.out=200)
    }
    nn <- length(theta.possible)
    objfn.check1 <- objfn.checkneg1 <- rep(0, nn)
     for(k in 1:nn) {
          objfn.check1[k] <- ProfileLogLikSQP(theta=theta.possible[k], gamma=1, nevents0=nevents0, nevents1=nevents1, 
                                             nrisk0=nrisk0, nrisk1=nrisk1, utimes0=utimes0, utimes1=utimes1,H0=H0, H1=H1)
          objfn.checkneg1[k] <- ProfileLogLikSQP(theta=theta.possible[k], gamma=-1, nevents0=nevents0, nevents1=nevents1, 
                                              nrisk0=nrisk0, nrisk1=nrisk1, utimes0=utimes0, utimes1=utimes1,H0=H0, H1=H1)
          print(k)
     }
     
     opt.gam1 <- min(objfn.check1)
     opt.gamneg1 <- min(objfn.checkneg1)
     opt.theta1 <- theta.possible[which.min(objfn.check1)]
     opt.thetaneg1 <- theta.possible[which.min(objfn.checkneg1)]
     best.gamma <- ifelse(opt.gam1 < opt.gamneg1, 1, -1)
     best.theta <- ifelse(best.gamma==1, opt.theta1, opt.thetaneg1)
     
   
  } else if(is.null(theta.fixed) & !is.null(gamma)) {
      theta.interval <- c(min(c(utimes0, utimes1)), .95*max(c(utimes0, utimes1)) )
      ## probably need to give more thought to the best choice for theta.interval
      if(gamma==1) {
           opt.gam1 <- optimize(f=ProfileLogLikSQP, interval=theta.interval, gamma=1, nevents0=nevents0, nevents1=nevents1, 
                                nrisk0=nrisk0, nrisk1=nrisk1, utimes0=utimes0, utimes1=utimes1,
                                H0=H0, H1=H1)
           theta.star <- opt.gam1$minimum
      } else if(gamma== - 1) {
           opt.gamneg1 <- optimize(f=ProfileLogLikSQP, interval=theta.interval, gamma=-1, nevents0=nevents0, nevents1=nevents1, 
                                  nrisk0=nrisk0, nrisk1=nrisk1, utimes0=utimes0, utimes1=utimes1,
                                  H0=H0, H1=H1)
           theta.star <- opt.gamneg1$minimum
      }
      best.gamma <- gamma
      best.theta <- theta.star
  } else if(!is.null(theta.fixed) & is.null(gamma)) {
      theta.possible <- NULL
      obj1 <- ProfileLogLikSQP(theta=theta.fixed, gamma = 1, nevents0, nevents1, nrisk0, nrisk1, utimes0, utimes1, H0, H1) 
      objneg1 <- ProfileLogLikSQP(theta=theta.fixed, gamma = -1, nevents0, nevents1, nrisk0, nrisk1, utimes0, utimes1, H0, H1) 
      
      best.gamma <- ifelse(obj1 < objneg1, 1, -1)
      best.theta <- theta.fixed
  } else if(!is.null(theta.fixed) & !is.null(gamma)) {
      theta.possible <- NULL
      best.theta <- theta.fixed
      best.gamma <- gamma
      objfn.check1 <- NULL
  }
  #cat("best.theta", best.theta, "\n")
  
  ## (4) 
  ## Change this function so that it can take a fixed value of gamma as well.
  tmp <- CumHazKnownTheta(best.theta, gamma=best.gamma, nevents0, nevents1, nrisk0, nrisk1, utimes0, utimes1,
                          H0=H0, H1=H1)

  if(is.null(objfn.check1)) {
     objfn.check1 <- objfn.checkneg1 <- tmp$DistNA
  }  
  hazard0 <- tmp$hazard0
  hazard1 <- tmp$hazard1
  surv0 <- cumprod(1 - hazard0)
  surv1 <- cumprod(1 - hazard1)
  ## note that this answer excludes the last jump point
  ans <- list(times0=tmp$utimes0, nevents0=tmp$nevents0, nrisk0=tmp$nrisk0, times1=tmp$utimes1, nevents1=tmp$nevents1,
              nrisk1=tmp$nrisk1, hazard0=hazard0, hazard1=hazard1, surv0=surv0, surv1=surv1, theta=best.theta,
              gamma=best.gamma, theta.possible=theta.possible, objfn.check1=objfn.check1, objfn.checkneg1=objfn.checkneg1,
              H0=H0, H1=H1, DistNA=tmp$DistNA)
  class(ans) <- "surv.delay"
  return(ans)
}