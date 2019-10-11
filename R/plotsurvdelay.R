plot.surv.delay <- function(obj, type="surv") {
  ## plotting function
  if(type=="surv") {
    plot(obj$times0, obj$surv0, type="n", las=1, xlab="Time", ylab="Survival Prob.",
         xlim=c(0, max(c(obj$times0, obj$times1))), ylim=c(0,1))
    lines(c(0,obj$times0), c(1, obj$surv0), type="s", lwd=2)
    lines(c(0,obj$times1), c(1, obj$surv1), type="s", col="red", lwd=2)
    abline(v=obj$theta, lwd=2)
    legend("bottomleft", legend=c("Control Arm", "Active Trt. Arm"), col=c("black", "red"), lwd=3,
           bty='n')
  } else if(type=="cum.hazard") {
    plot(obj$times0, cumsum(obj$hazard0), type="n", las=1, xlab="Time", ylab="Cumulative Hazard",
         xlim=c(0, max(c(obj$times0, obj$times1))))
    lines(obj$times0, cumsum(obj$hazard0), type="s", lwd=2)
    lines(obj$times1, cumsum(obj$hazard1), type="s", col="red", lwd=2)
    abline(v=obj$theta, lwd=2)
    legend("topleft", legend=c("Control Arm", "Active Trt. Arm"), col=c("black", "red"), lwd=3,
           bty='n')
  }
}