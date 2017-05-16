\name{excessm}
\alias{excessm}
\alias{localmax}
\alias{lambdaweight}
\title{Excess Mass Function
}
\description{Algorithm which calculates the empirical excess mass for a given \eqn{\lambda} and given maximal number of modes.
}
\usage{
excessm(x, lambda, M = 1, UpToM = FALSE)
}
\arguments{
  \item{x}{data in form of a vector 
}
  \item{lambda}{\eqn{\lambda} of interest 
}
  \item{M}{maximal number of modes 
}
  \item{UpToM}{if true, the intervals for modes up to M are returned 
}
}

\value{
\item{intervals}{Matrix containing the empirical \eqn{\lambda}-clusters. Line \eqn{i} consists of information about the \eqn{i.-\lambda}-cluster, where \eqn{i \le M}. First (Second) column contains the corresponding starting (ending) sorted observation index. Third (Fourth) column the starting (ending) value. In case \eqn{UpToM} is true and \eqn{M>1}, a list is returned with the intervals for all modes \eqn{\le M}}
\item{excess_mass}{returns a vector with excess masses, the \eqn{i.} entry is the excess mass achieved with \eqn{i} modes}
}
\references{Muller, D. W. and Sawitzki, G., 09.1991, "Excess Mass Estimates and Tests for Multimodality", \emph{Journal of the American Statistical Association} \bold{, Vol. 86, No. 415, pp. 738--746}, \url{http://www.jstor.org/stable/2290406}
}
\author{Marc-Daniel Mildenberger \email{mildenberger.stat@web.de}, based on earlier code from Dr. Guenther Sawitzki \email{gs@statlab.uni-heidelberg.de}
}
\note{Please note that an allowance for \eqn{M} modes does not necessarily result in \eqn{M} \eqn{\lambda}-clusters. Hence, the number of intervals returned can be smaller than \eqn{M}. In this case a warning will be displayed. The vector \eqn{excess_mass} does have less than \eqn{M} entries.
}

\seealso{\code{\link{exmplot}, \link{exmsilhouette}, \link{mexmsilhouette}}
}
\examples{
library(MASS)
attach(geyser)

##calculating excess mass for duration of 'Old Faithful Geyser' for lambda=0.2 allowing for one mode
excessm(duration, lambda=0.2)

##same as above, but allowing for up to three modes
excessm(duration, lambda=0.2, M=3) 

#returns the intervals for modes 1,2 and 3
excessm(duration, lambda=0.2, M=3, UpToM=TRUE) 
}
\keyword{ htest }
\keyword{ univ }
\keyword{ nonparametric }