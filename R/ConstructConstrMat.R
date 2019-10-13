ConstructConstrMat <- function(times0, times1, theta) {
  ### Function for constructing the main constraint matrix C
  #   we should have Cw \geq 0
  #   This is the constraint matrix under the assumption that
  #   the control arm initially "dominates" the active treatment arms 
  
  times.tot <- sort(c(times0, times1))
  
  
  tau <- length(times.tot)
  A <- matrix(0, nrow=tau, ncol=length(times0))
  B <- matrix(0, nrow=tau, ncol=length(times1))
  Kstar <- max(which(times.tot < theta))
  for(j in 1:tau) {
    if(j <= Kstar) {
      A[j,] <- ifelse(times0 <= times.tot[j], -1, 0)
      B[j,] <- ifelse(times1 <= times.tot[j], 1, 0)
    } else {
      A[j,] <- ifelse(times0 <= times.tot[j], 1, 0)
      B[j,] <- ifelse(times1 <= times.tot[j], -1, 0)
    }
  }
  ans <- cbind(A, B)
  return(ans)
}