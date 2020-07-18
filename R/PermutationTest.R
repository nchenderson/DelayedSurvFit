PermutationTest <- function(times, events, trt, nperms=100, type="KM") {
  
    n <- length(trt)
    WW <- rep(0, nperms)
    pi0 <- mean(trt==0)
    ds.obj <- DelayedSurvFit(times, events, trt) 
    W.obs <- WeightedKMStat(ds.obj$theta, ds.obj$discretized.times, events, trt) 
    for(k in 1:nperms) {
       trt.assign <- sample(0:1, size=n, replace=TRUE, prob=c(pi0, 1 - pi0))
       ds.obj <- DelayedSurvFit(times, events, trt.assign, final.iter=5) ## change inner.iter and final.iter here (at least final.iter=5) (maybe inner.iter=20)
       ## Compute weighted Kaplan-Meier statistic 
       WW[k] <- WeightedKMStat(ds.obj$theta, ds.obj$discretized.times, events, trt.assign) 
    }
    perm.pval <- (1 + sum(WW >= W.obs))/(nperms + 1)
    print(WW)
    print(W.obs)
    return(perm.pval)
}