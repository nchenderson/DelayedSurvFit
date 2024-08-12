\name{DelayedSurvFit}
\alias{DelayedSurvFit}
\title{
Survival functions with single-crossing constraints
}
\description{
This computes nonparametric estimates of the survival curves
in two treatment arms under the constraint that the two survival curves can cross
at most one time.
}
\usage{
DelayedSurvFit(times, events, trt, gamma.fixed=NULL, theta.fixed=NULL, 
               weights=NULL, max.times=100, inner.iter=50, final.iter=1000, 
               verbose=TRUE)
}
\arguments{
  \item{times}{
vector of follow-up times
}
\item{events}{
vector of event indicators: 1 for observed events, 0 for censored observations
}
\item{trt}{
vector of treatment arm assignment: 1 for active treatment arm, 0 for control arm
}
\item{gamma.fixed}{
optional fixed value of the superior long-term survival parameter: gamma.fixed = 1 means the active treatment arm has better long-term survival while gamma.fixed = -1 means the control arm has better long-term survival. This is estimated if not provided by the user. 
}
\item{theta.fixed}{
optional fixed value of theta. This is estimated if not provided by the user.
}
\item{weights}{
optional vector of individual-level weights. This should be a vector with the same length as \code{times}.
}
\item{max.times}{
The maximum number of bins to use for discretization of the follow-up times.
}
\item{inner.iter}{
The number of iterations used to find the crossing time parameters.
}
\item{final.iter}{
Given estimated crossing time parameter estimates, number of iterations used to find survival probabilities.
}
\item{verbose}{
If \code{TRUE}, progress of the algorithm is printed to the screen.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
an object of class \code{surv.delay}

This is a list withthe following components
  \item{times}{Vector of length m which has the unique, sorted follow-times after discretization.  }
  \item{surv0}{Vector of length m which contains estimated control arm survival probabilities at each of the time points in the vector \code{times}}
  \item{surv1}{Vector of length m which contains estimated active treatment arm survival probabilities at each of the time points in the vector \code{times}}
  \item{haz0}{Vector of length m which contains estimated control arm hazards at each of the time points in the vector \code{times}}
  \item{haz1}{Vector of length m which contains estimated active treatment arm hazards at each of the time points in the vector \code{times}}
  \item{nevents0}{Vector of length m which contains estimated control arm survival probabilities at each of the time points in the vector \code{times}}
  \item{nrisk0}{Vector of the number at risk in control arm}
  \item{nevents1}{Vector of number of events in active treatment arm}
  \item{nrisk1}{Vector of the number at risk in active treatment arm}
  \item{theta}{Estimated value of the survival curve crossing time}
  \item{gamma}{Estimated value of which arm has superior long-term survival}
  \item{discretized.times}{Vector of discretized follow-up times used in the computation}
  \item{negloglik.val}{Value of the negative log-likelihood at convergence}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Nicholas Henderson
}
\note{
%%  ~~further notes~~
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(5172)

### Generate simulated survival data
n <- 400
cens.time <- 8
X1 <- rweibull(n, shape=2.3, scale=7)
X2 <- rweibull(n, shape=.8, scale=10)
Y1 <- pmin(X1, cens.time)
Y2 <- pmin(X2, cens.time)
e1 <- as.numeric(X1 < cens.time)
e2 <- as.numeric(X2 < cens.time)

times <- c(Y1, Y2)  ## follow-up times
events <- c(e1, e2) ## event indicators
trt <- rep(c(0,1), each=n)  ## treatment arm assignment

## Run the single-crossing constrained estimation procedure with a fixed 
## value of gamma = -1 and theta = 5.8
dsf.obj <- DelayedSurvFit(times, events, trt, gamma.fixed=-1, 
                          theta.fixed=5.8, max.times=100) 

km0 <- cumprod(1 - dsf.obj$nevents0/dsf.obj$nrisk0)
km1 <- cumprod(1 - dsf.obj$nevents1/dsf.obj$nrisk1)

## Plot the single-crossing constrained estimate with the KM estimate:
plot(dsf.obj, xlab="time", ylab="survival probability")
lines(dsf.obj$times, km0, type="s", lty=2)
lines(dsf.obj$times, km1, type="s", lty=2, col="red")


}
\keyword{ survival }
\keyword{ nonparametric }
