DiscretizeTimes <- function(times, max.times=max.times) {
   chc <- hclust(dist(times))
   memb <- cutree(chc, k=max.times)
   new.times <- rep(0, length(times))
   for(k in 1:100) {
     new.times[memb==k] <- median(times[memb==k])
   }
   return(new.times)
}