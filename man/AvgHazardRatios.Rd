\name{AvgHazardRatios}
\alias{AvgHazardRatios}
%- Also NEED an '\alias' for EACH other topic documented here.
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
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
a <- 2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
\keyword{ nonparametric }
