rrml <- function(obj, t0, t1, tau0=NULL, tau1=NULL)
{
  if(is.null(tau0)) {
    tau0 <- max(obj$times)
  }
  if(is.null(tau1)) {
    tau1 <- max(obj$times)
  }
  ss1 <- SurvFn(obj, arm=1)
  ss0 <- SurvFn(obj, arm=0)
  rmst.vals.tau <- rmst(obj, tau0=tau0, tau1=tau1)
  rmst.vals.t <- rmst(obj, tau0=t0, tau1=t1)
  
  rmst.diff0 <- rmst.vals.tau$rmst0 - rmst.vals.t$rmst0
  rmst.diff1 <- rmst.vals.tau$rmst1 - rmst.vals.t$rmst1
  rrml0 <- rmst.diff0/ss0(t0)
  rrml1 <- rmst.diff1/ss1(t1)
 
  return(list(rrml0=rrml0, rrml1=rrml1))
}