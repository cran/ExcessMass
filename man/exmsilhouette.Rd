\name{exmsilhouette}
\alias{exmsilhouette}
\alias{colorizing}
\title{(Single) Excess Mass Silhouette
}
\description{Produces an excess mass plot and the corresponding numerical values if required. 
}
\usage{
exmsilhouette(xdata, M = 1, CutOff = 1,steps = 50,rug = TRUE,
Lambda = NULL,col = FALSE,rdata = FALSE,label = TRUE)
}
\arguments{
  \item{xdata}{data in form of vector
}
  \item{M}{the maximal number of modes
}
  \item{CutOff}{determines the cut off value and hence the appearance of the graph
}
  \item{steps}{number of different \eqn{\lambda}s which are calculated / drawn
}
  \item{rug}{draws a rug plot at the bottom of the graph
}
  \item{Lambda}{allows to specify an own vector of \eqn{\lambda}s which are drawn in the graph
}
  \item{col}{lines get colored in purple (\eqn{em > 0.75}), green (\eqn{0.75 \ge em > 0.5}), turquoise (\eqn{0.5 \ge em > 0.25}), blue (\eqn{0.25 \ge em > 0.05}) and black (\eqn{0.05 \ge em})
}
  \item{rdata}{a numerical output is returned
}
  \item{label}{allows to reduce labeling
}
}
\details{\eqn{CutOff} should not be set too small or too large, as this results in meaningless graphs. 
The excess mass for several \eqn{\lambda}s can be calculated by specifying the Lambda. 
}
\value{A plot is always produced. By setting \eqn{rdata=TRUE} numerical results are returned in form of a two-dimensional list. The first argument specifies \eqn{\lambda}. This means that if \eqn{steps=50} with \eqn{1} (\eqn{50}) you get access to the numerical results for the smallest (largest) \eqn{\lambda}.

In case no \eqn{Lambda} vector is used for each \eqn{\lambda}, the following information is displayed:
\item{\eqn{[,1]}}{value of \eqn{\lambda},}
\item{\eqn{[,2]}}{ calculated \eqn{\lambda}-clusters}
\item{\eqn{[,3]}}{excess mass vector}
The last two components are presented in the way known from  the excess mass function. In case \eqn{Lambda} was set manually the value of \eqn{\lambda} is not returned, as it is known.
}
\references{Muller, D. W. and Sawitzki, G., 09.1991, "Excess Mass Estimates and Tests for Multimodality", \emph{Journal of the American Statistical Association} \bold{, Vol. 86, No. 415, pp. 738--746}, \url{http://www.jstor.org/stable/2290406}
}
\author{Marc-Daniel Mildenberger \email{mildenberger.stat@web.de}, based on earlier code from Dr. Guenther Sawitzki \email{gs@statlab.uni-heidelberg.de} 
}

\seealso{\code{\link{excessm}, \link{mexmsilhouette}, \link{exmplot}}
}
\examples{
library(MASS)
attach(geyser)

##Plot allowing for up to two modes and reduced CutOff value
exmsilhouette(duration, M=2, CutOff=1.25)

##Plot with twice the default number of steps, omitting rug plot, 
##colorizing the graph and asking for numerical output
res <- exmsilhouette(duration, M=2, CutOff=1.25, steps=100, rug=FALSE, col=TRUE, rdata=TRUE)

##Specifying Lambda and requesting numerical output
L=seq(.01,0.25,0.005)
res <- exmsilhouette(duration, M=3, Lambda=L, col=TRUE, rdata=TRUE)
}

\keyword{ hplot }
\keyword{ htest }
\keyword{ univ }
\keyword{ nonparametric }
