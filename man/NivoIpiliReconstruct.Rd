\name{NivoIpiliReconstruct}
\alias{NivoIpiliReconstruct}
\docType{data}
\title{
Approximate survival data for Nivolumab plus Ipilimumab
}
\description{
Reconstructed survival data from a study involving
the comparison of Nivolumab plus Ipilimumab vs. chemotherapy
in the treatment of patients having advanced non–small-cell lung cancer.
}
\usage{data("NivoIpiliReconstruct")}
\format{
  A data frame with 1166 observations on the following 3 variables.
  \describe{
    \item{\code{times}}{vector of follow-up times}
    \item{\code{evnt}}{vector of event indicators}
    \item{\code{trt}}{treatment assignment indicators; 1 - Nivolumab + Ipili, 0 - chemotherapy}
  }
}
\details{
Approximate survival data reconstructed from survival curves extracted
from Hellmann et al. (2019)
}
\source{
Hellmann, Matthew D., et al. (2019). "Nivolumab plus Ipilimumab in Advanced Non–Small-Cell Lung Cancer." New England Journal of Medicine.
}
\references{
Hellmann, Matthew D., et al. (2019). "Nivolumab plus Ipilimumab in Advanced Non–Small-Cell Lung Cancer." New England Journal of Medicine.
}
\examples{
data(NivoIpiliReconstruct)
head(NivoIpili) 
}
\keyword{datasets}
