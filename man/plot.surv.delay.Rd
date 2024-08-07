\name{plot.surv.delay}
\alias{plot.surv.delay}
\title{
Plot method for \code{surv.delay} objects
}
\description{
A plot of either the estimated survival functions or hazard functions is produced.
If the "surv" type is selected, the survival functions for each treatment arm is produced.
If the "hazard" type is selected, the hazard functions for each treatment arm is produced.
}
\usage{
plot.surv.delay(x, ..., type="surv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{surv.fit}, usually returned by either the \code{DelayedSurvFit} or \code{DelayedHazFit} functions
}
  \item{type}{
The type of plot desired. Should be "surv" if the survival functions are desired and should be "hazard" 
if the hazard functions are desired.
}
\item{\dots}{
other arguments that will be passed forward to the default plot method such as \code{ylim} or \code{xlab}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A \code{data.frame} with the following variables
   \item{x}{the values of the x coordinates used in the plot}
   \item{y}{the values of the y coordinates used in the plot}
   \item{trt}{an indicator of which treatment arm the (x,y) coordinate comes from}
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
\dontrun{
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

## Plot fitted survival functions
plot(dfs.obj, xlab="Survival Prob.", ylab="Time", las=1)

## Plot fitted hazard functions
plot(dfs.obj, type="hazard", xlab="Hazard", ylab="Time", las=1)
}
}
\keyword{ aplot }
\keyword{ survival }
