\name{DelayedSurvFit}
\alias{DelayedSurvFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting Survival Curves with Potential Crossings
}
\description{
aa
}
\usage{
DelayedSurvFit(times, events, trt, gamma, theta.fixed=NULL)
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
a = 2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
