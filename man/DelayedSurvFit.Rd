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
DelayedSurvFit(times, events, trt, gamma=NULL, theta.fixed=NULL, max.times=100,
               inner.iter=50, final.iter=1000, verbose=TRUE)
}
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
an object of class surv.delay
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

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

dsf.obj <- DelayedSurvFit(times, events, trt, gamma=-1, theta.fixed=5.8, max.times=100) 

km0 <- cumprod(1 - dsf.obj$nevents0/dsf.obj$nrisk0)
km1 <- cumprod(1 - dsf.obj$nevents1/dsf.obj$nrisk1)

plot(dsf.obj)
lines(dsf.obj$times, km0, type="s", lty=2)
lines(dsf.obj$times, km1, type="s", lty=2, col="red")


}
\keyword{ survival }
\keyword{ nonparametric }
