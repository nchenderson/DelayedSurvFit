plot.surv.delay <- function(x, ..., type="surv") {
  ## plotting function
  if(type=="surv") {
     plot(c(0, rep(x$times, 2)), c(1,x$surv0, x$surv1), type="n", ...)
     lines(c(0, x$times), c(1, x$surv0), type="s", lwd=2)
     lines(c(0, x$times), c(1, x$surv1), type="s", col="red", lwd=2)
     abline(v=x$theta, lwd=2)
     legend("bottomleft", legend=c("Control Arm", "Active Trt. Arm"), col=c("black", "red"), lwd=3,
           bty='n')
     ans <- data.frame(x=rep(x$times, 2), y=c(x$surv0, x$surv1), trt=rep(c(0,1), each=length(x$times)))
  } else if(type=="hazard") {
     plot(rep(x$times, 2), c(x$haz0, x$haz1), type="n", ...)
     lines(x$times, x$haz0, type="s", lwd=2)
     lines(x$times, x$haz1, type="s", col="red", lwd=2)
     abline(v=x$theta, lwd=2)
     legend("topleft", legend=c("Control Arm", "Active Trt. Arm"), col=c("black", "red"), lwd=3,
            bty='n')
     ans <- data.frame(x=rep(x$times, 2), y=c(x$haz0, x$haz1), trt=rep(c(0,1), each=length(x$times)))
  }
  return(ans)
}
