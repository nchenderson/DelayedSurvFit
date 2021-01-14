AvgHazardRatios <- function(obj) {
  
   if(obj$theta > min(obj$times) && obj$theta <= max(obj$times[-length(obj$times)])) {
      lt.theta <- obj$times <= obj$theta
      time.delta <- diff(c(0, obj$times))
      haz.ratios <- obj$haz1/(obj$haz0 + obj$haz1)  
      pre.avg <- sum(haz.ratios*time.delta*lt.theta)/obj$theta
      post.avg <- sum(haz.ratios*time.delta*(1 - lt.theta))/(max(obj$times) - obj$theta)
   } else {
      time.delta <- diff(obj$times)
      haz.ratios <- obj$haz1[-length(obj$times)]/(obj$haz0[-length(obj$times)] + obj$haz1[-length(obj$times)])  
      pre.avg <- sum(haz.ratios*time.delta)
      post.avg <- NULL
   }
   return(list(precross.avg=pre.avg, postcross.avg=post.avg))
}


