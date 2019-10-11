ConstructConstrMat <- function(times0, times1, Kstar) {
  tau <- length(times1)
  A <- matrix(0, nrow=tau, ncol=length(times0))
  B <- matrix(0, nrow=tau, ncol=tau)
  for(j in 1:tau) {
    if(j <= Kstar) {
      A[j,] <- ifelse(times0 <= times1[j], -1, 0)
      B[j,] <- ifelse(times1 <= times1[j], 1, 0)
    } else {
      A[j,] <- ifelse(times0 <= times1[j], 1, 0)
      B[j,] <- ifelse(times1 <= times1[j], -1, 0)
    }
  }
  ans <- cbind(A, B)
  return(ans)
}