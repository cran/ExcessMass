\name{exmplot}
\alias{exmplot}
\title{Excess Mass Lambda Plot
}
\description{Produces an excess mass lambda plot and calculates the maximal excess mass difference achieved by allowing for an additional mode.
}
\usage{
exmplot(xdata, M=1, CutOff=1, steps=50, Lambda=NULL)
}
\arguments{
  \item{xdata}{data in form of a vector
}
  \item{M}{the maximal number of modes
}
  \item{CutOff}{determines the cut off value and hence the level up to which the \eqn{\lambda}s are calculated
}
  \item{steps}{number of different \eqn{\lambda}s which are calculated / drawn
}
  \item{Lambda}{allows specifying an own vector of \eqn{\lambda}s
}
}
\details{\eqn{CutOff} should not be set too small or too large, as this results in meaningless graphs. 
The excess mass for several \eqn{\lambda}s can be calculated by specifying the Lambda. 
}
\value{An Excess Mass Lambda plot is produced. The lines in the plot are sorted by the maximal number of modes from left to right, due to the monotonicity of the excess mass in M.
\item{max_dist}{The i. entry is the maximal distance of the excess mass by allowing for up to i+1 instead of i modes}
\item{max_dist_Lambda}{Shows the \eqn{\lambda} at which the maximal excess mass difference is achieved}
}
\references{Muller, D. W. and Sawitzki, G., 09.1991, "Excess Mass Estimates and Tests for Multimodality", \emph{Journal of the American Statistical Association} \bold{, Vol. 86, No. 415, pp. 738--746}, \url{http://www.jstor.org/stable/2290406}
}
\author{Marc-Daniel Mildenberger \email{mildenberger.stat@web.de}, based on earlier code from Dr. Guenther Sawitzki \email{gs@statlab.uni-heidelberg.de} 
}

\seealso{\code{\link{excessm}, \link{exmsilhouette}, \link{mexmsilhouette}}
}
\examples{
library(MASS)
attach(geyser)

##calculating the maximal excess mass difference for duration of 'Old Faithful Geyser' for M=3
exmplot(duration, M=3)

##Plotting the excess mass against lambda for modes 1-5, 
##increase CutOff value, double the number of steps
exmplot(duration, M=5, CutOff=1.2, steps=100)

##Specifying Lambda
Lambda=seq(.0,0.5,0.005)
exmplot(duration, M=7, Lambda=Lambda)
}

\keyword{ hplot }
\keyword{ htest }
\keyword{ univ }
\keyword{ nonparametric }
