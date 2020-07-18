

WeightedKMStat <- function(cross.point, times, events, trt) {
    sfit.cens <- survfit(Surv(times, 1 - events) ~ trt)
    sfit <- survfit(Surv(times, events) ~ trt)
    risk.data.cens <- data.frame(strata = summary(sfit.cens, times = times, extend = TRUE)$strata,
                            time = summary(sfit.cens, times = times, extend = TRUE)$time,
                            surv = summary(sfit.cens, times=times, extend=TRUE)$surv)
    risk.data <- data.frame(strata = summary(sfit, times = times, extend = TRUE)$strata,
                                 time = summary(sfit, times = times, extend = TRUE)$time,
                                 surv = summary(sfit, times=times, extend=TRUE)$surv)
    
    S1 <- risk.data$surv[risk.data$strata=="trt=1"]
    S0 <- risk.data$surv[risk.data$strata=="trt=0"]
    tt <- risk.data$time[risk.data$strata=="trt=0"]
    
    G1 <- risk.data.cens$surv[risk.data.cens$strata=="trt=1"]
    G0 <- risk.data.cens$surv[risk.data.cens$strata=="trt=0"]
    tt <- risk.data.cens$time[risk.data.cens$strata=="trt=0"]
    
    n0 <- sum(trt==0)
    n1 <- sum(trt==1)
    nn <- n0 + n1
    
    #Gfn <- stepfun(tt, c(1, (nn*G1*G0)/(n0*G0 + n1*G1))*(c(0,tt) >= cross.point)*c(0, S1 - S0)) 
    #plot(Gfn)
    FF <- c(1, (nn*G1*G0)/(n0*G0 + n1*G1))*(c(0,tt) >= cross.point)*c(0, S1 - S0)
    FF <- FF[!is.na(FF)]
    tau <- length(FF)
    W.stat <- sum(c(tt[1], diff(tt[1:(tau-1)]))*FF[-tau])
    #integrate(Gfn, lower=0, upper=8, subdivisions=500L)
    return(W.stat)
}


