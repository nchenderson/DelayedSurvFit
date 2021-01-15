\name{SurvFn}
\alias{SurvFn}
\title{
Arm-specific survival functions
}
\description{
Given a fitted object from either \code{DelayedSurvFit} or \code{DelayedHazFit}, this
will return a step function which represents the estimated survival function
for either the control arm or active treatment arm.
}
\usage{
SurvFn(obj, arm = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
An object of class \code{surv.delay}
}
  \item{arm}{
The treatment arm; \code{arm = 0} for the control arm; \code{arm = 1} for the active treatment arm
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A function of class "stepfun"
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

dfs.obj <- DelayedSurvFit(times, events, trt, max.times=50) 

S1 <- SurvFn(dfs.obj, arm=1) # active treatment arm survival function
S0 <- SurvFn(dfs.obj, arm=0) # control treatment arm survival function

## Difference in survival function at time 2
S1(2) - S0(2)

## Difference in survival function at time 5
S1(5) - S0(5)
}
\keyword{ survival }
\keyword{ nonparametric }
