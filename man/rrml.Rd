\name{rrml}
\alias{rrml}
\title{
Restricted residual mean life
}
\description{
Computes the restricted residual mean life for both treatment arms.
}
\usage{
rrml(obj, t0, t1, tau0 = NULL, tau1 = NULL)
}
\arguments{
  \item{obj}{
An object of class \code{surv.delay}
}
  \item{t0}{
The lower trunctation time in the control arm.
}
  \item{t1}{
The lower trunctation time in the active treatment arm.
}
  \item{tau0}{
The upper trunctation time in the control arm.
}
  \item{tau1}{
The upper trunctation time in the active treatment arm.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with the following two components
  \item{rrml0}{The restricted residual mean life in the control arm}
  \item{rrml1}{The restricted residual mean life in the active treatment arm}
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
rrml(dfs.obj, t0=3, t1=3)
}
}
\keyword{ survival }
\keyword{ robust }
