\name{searchMaxLambda}
\alias{searchMaxLambda}
\title{Searching For A Maximal Lambda
}
\description{Function which gives a rough approximation of maximal \eqn{\lambda}.
}
\usage{
searchMaxLambda(x, limcount = 4, step = 1.05, trylambda = 0.01)
}
\arguments{
  \item{x}{data in form of a vector
}
  \item{limcount}{divided by the square root of the number of data points. The result determines the cut off value.
}
  \item{step}{determines step size
}
  \item{trylambda}{initial \eqn{\lambda}
}
}
\details{Excess mass is calculated for \eqn{trylambda}. In case the resulting excess mass is larger (smaller) than the cut off value, \eqn{trylambda} is set as \eqn{trylambda * step} (respectively \eqn{trylambda / step}) and excess mass is calculated again until it is smaller (larger) than the cut off value. The corresponding \eqn{trylambda} is returned.


The approximation is done allowing only for one \eqn{\lambda}-cluster, as scans including more \eqn{\lambda}-clusters have high computational costs due to the recursive structure of the algorithm. 
}
\value{\eqn{trylambda} calculated as described in Details.
}
\references{Muller, D. W. and Sawitzki, G., 09.1991, "Excess Mass Estimates and Tests for Multimodality", \emph{Journal of the American Statistical Association} \bold{, Vol. 86, No. 415, pp. 738--746}, \url{http://www.jstor.org/stable/2290406}
}
\author{Marc-Daniel Mildenberger \email{mildenberger.stat@web.de}, based on earlier code from Dr. Guenther Sawitzki \email{gs@statlab.uni-heidelberg.de} 
}

\seealso{\code{\link{excessm}, \link{exmplot}, \link{exmsilhouette}, \link{mexmsilhouette}
}}
\examples{
library(MASS)
attach(geyser)

#Calculating Lambda using standard settings
searchMaxLambda(duration)

#Calculating Lambda, reducing cut off value and step. Setting another initial lambda
searchMaxLambda(duration, limcount = 5, step = 1.01, trylambda = 1)
}

\keyword{ htest }
