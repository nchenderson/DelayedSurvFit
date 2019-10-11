DelayedHazard <- function(obj, bw=NULL) {
  ## bw - bandwidth
  if(is.null(bw)) {
    bw0 <- bw.nrd(obj$times0)
    bw1 <- bw.nrd(obj$times1)
    ## need to add a better bandwidth selection later
  } else {
    bw0 <- bw1 <- bw
  }
  grid.points0 <- seq(min(obj$times0), max(obj$times0), length.out=100)
  grid.points1 <- seq(min(obj$times1), max(obj$times1), length.out=100)
  
  Kmat0 <- dnorm(outer(grid.points0, obj$times0, FUN="-"), sd=bw0)
  Kmat1 <- dnorm(outer(grid.points1, obj$times1, FUN="-"), sd=bw1)
  
  smooth.hazard0 <- as.vector(Kmat0%*%obj$hazard0)
  smooth.hazard1 <- as.vector(Kmat1%*%obj$hazard1)
  
  ymax <- max(c(smooth.hazard0, smooth.hazard1))
  xmax <- max(max(obj$times0), max(obj$times1))
  plot(grid.points0, smooth.hazard0, type="n", las=1, xlab="Time", ylab="hazard", 
       main="Estimated hazard functions", ylim=c(0, ymax), xlim=c(0, xmax))
  lines(grid.points0, smooth.hazard0, lwd=2)
  lines(grid.points1, smooth.hazard1, col="red", lwd=2)
  legend("topleft", legend=c("Control Arm", "Active Trt. Arm"), col=c("black", "red"), lwd=3,
         bty='n')
  
  return(list(time0=grid.points0, time1=grid.points1, smooth.hazard0=smooth.hazard0,
              smooth.hazard1=smooth.hazard1))
}