FindInitialVectors <- function(theta, utimes0, utimes1) {
  eps <- 0.01
  w0 <- rep(0, length(utimes0))
  w1 <- rep(0, length(utimes1))
  w1[utimes1 < theta] <- rep(eps, sum(utimes1 < theta))
  n00 <- sum(utimes0 < theta)
  w0[utimes0 < theta] <- rep(eps/n00, n00)
  
  gap <- sum(w1) - sum(w0)
  w0[utimes0 >= theta] <- rep(eps, sum(utimes0 >= theta))
  n11 <- sum(utimes1 >= theta)
  w1[utimes1 >= theta] <- rep(eps/n11, n11)
  f.pass <- min(which(utimes0>=theta))
  w0[f.pass] <- gap + eps
  ans <- list(w0=w0, w1=w1)
  return(ans)
}