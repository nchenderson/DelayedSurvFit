FindActiveSet <- function(theta, gamma, utimes0, utimes1) {
    ## This for the case when control treatment dominates active treatment early on
    ## In this case, we need to have w_{0j} = 0 for j such that t_{0j} \leq t_{11}
    ## We also want to ignore w_{1j} such that \theta < w_{1j} < \min\{t_{0j}: t_{0j} > theta}
  
    if(gamma==1) {
        ind0 <- utimes0 <= utimes1[1]
        ## w0[ind0] = 0 
        tmp.ind <- utimes0 > theta
        if(sum(tmp.ind) > 0) {
             first.pass.theta <- min(utimes0[utimes0 > theta]) # make this more robust later
        } else {
             first.pass.theta <- Inf
        }
        ind1 <- utimes1 > theta & utimes1 <= first.pass.theta
  
        ans <- list(active.set0 = !ind0, active.set1=!ind1)
    } else if(gamma == -1) {
        ind1 <- utimes1 <= utimes0[1]
        
        tmp.ind <- utimes1 > theta
        if(sum(tmp.ind) > 0) {
          first.pass.theta <- min(utimes1[utimes1 > theta]) # make this more robust later
        } else {
          first.pass.theta <- Inf
        }
       # first.pass.theta <- min(utimes1[utimes1 > theta])
        ind0 <- utimes0 > theta & utimes0 <= first.pass.theta
        
        ans <- list(active.set0 = !ind0, active.set1=!ind1)
    }
    return(ans)
}



