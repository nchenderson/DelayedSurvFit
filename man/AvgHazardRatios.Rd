\name{AvgHazardRatios}
\alias{AvgHazardRatios}
\title{
Pre and post-crossing average hazard ratios
}
\description{
This function computes the average treatment-to-total hazard ratio
for the time periods before and after a crossing occurs. 
}
\usage{
AvgHazardRatios(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
An object of class \code{surv.delay}
}
}
\details{
The crossing time can either represent a crossing of the survival functions
or a crossing a crossing of the hazard functions. If one wants to use
a crossing time of the survival functions, make sure the 
\code{surv.delay} object was obtained from running the \code{DelayedSurvFit}
function. Likewise, if one wants to use
a crossing time of the hazard functions, make sure the 
\code{surv.delay} object was obtained from running the \code{DelayedHazFit}
function. 
}
\value{
A list with the following two components
  \item{precross.avg}{Average treatment-to-total hazard ratio before the crossing time}
  \item{postcross.avg}{Average treatment-to-total hazard ratio after the crossing time}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Nicholas Henderson
}
\note{
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

dhz.obj <- DelayedHazFit(times, events, trt, max.times=50) 

AvgHazardRatios(dhz.obj)
}
\keyword{ survival }
\keyword{ nonparametric }
