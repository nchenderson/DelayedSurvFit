SurvFn <- function(obj, arm=1) {
    if(arm==1) {
        SS <- stepfun(obj$times, c(1, obj$surv1), right=TRUE)
    } else if (arm==0) {
        SS <- stepfun(obj$times, c(1, obj$surv0), right=TRUE)
    }
    return(SS)
}