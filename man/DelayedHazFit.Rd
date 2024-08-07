\name{DelayedHazFit}
\alias{DelayedHazFit}
\title{
Hazard functions with single-crossing constraints
}
\description{
This computes nonparametric estimates of discrete hazards 
from survival data in two treatment arms under
the constraint that the two hazard functions can cross
at most one time.
}
\usage{
DelayedHazFit(times, events, trt, gamma = NULL, theta.fixed = NULL, 
max.times = 100, inner.iter = 50, final.iter = 1000, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{
vector of follow-up times
}
\item{events}{
vector of event indicators
}
\item{trt}{
vector of treatment arm assignment
}
\item{gamma}{
value of the gamma parameter
}
\item{theta.fixed}{
optional fixed value of theta. Estimated if not provided.
}
\item{max.times}{
tbf
}
\item{inner.iter}{
tbf
}
\item{final.iter}{
tbf
}
\item{verbose}{
tbf
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
an object of class "surv.delay"
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

dhz.obj <- DelayedHazFit(times, events, trt, max.times=50) 

plot(dhz.obj$times, dhz.obj$haz0, type="n")
ntimes <- length(dhz.obj$times)
for(k in 1:ntimes) {
  if(dhz.obj$haz0[k] > dhz.obj$haz1[k]) {
     lines(c(dhz.obj$times[k], dhz.obj$times[k]), c(0, dhz.obj$haz0[k]), lwd=2)
     lines(c(dhz.obj$times[k], dhz.obj$times[k]), c(0, dhz.obj$haz1[k]), col="red", lwd=2)
  } else {
     lines(c(dhz.obj$times[k], dhz.obj$times[k]), c(0, dhz.obj$haz1[k]), col="red", lwd=2)
     lines(c(dhz.obj$times[k], dhz.obj$times[k]), c(0, dhz.obj$haz0[k]), lwd=2)
  }
}
abline(v=dhz.obj$theta, lwd=2, lty=2)
}
\keyword{ survival }
\keyword{ nonparametric }
