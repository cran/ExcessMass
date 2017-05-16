\name{mexmsilhouette}
\alias{mexmsilhouette}
\title{(Multiple) Excess Mass Silhouette
}
\description{Produces a graph with several excess mass plots allowing for different maximal numbers of modes/ cut off values.
}
\usage{
mexmsilhouette(xdata, M = 1:3, CutOff = c(1,2,5), steps = 30,
Lambda = NULL, col = FALSE, rug = TRUE, rdata = FALSE)
}
\arguments{
  \item{xdata}{data in form of a vector
}
  \item{M}{vector containing the max. number of modes
}
  \item{CutOff}{vector which determines the cut off values and hence the appearance of the graph
}
  \item{steps}{number of different \eqn{\lambda}s which are calculated / drawn
}
  \item{Lambda}{allows to specify an own vector of \eqn{\lambda}s which is drawn in the graph
}
  \item{col}{lines get colored in purple (\eqn{em > 0.75}), green (\eqn{0.75 \ge em > 0.5}), turquoise (\eqn{0.5 \ge em > 0.25}), blue (\eqn{0.25 \ge em > 0.05}) and black (\eqn{0.05 \ge em})
}
  \item{rug}{draws a rug plot at the bottom of the graph
}
  \item{rdata}{a numerical output is returned
}
}
\details{\eqn{CutOff} should not be set too small or too large, as this results in meaningless graphs.
}
\value{Always a graph with multiple plots is produced. Each column contains another maximal number of modes and each row another CutOff factor.

Setting rdata=TRUE numerical results are returned in form of a list. If the number of modes and the CutOff parameter contain just one element, the output of "mexmsilhouette" and "exmplot" are equal. 

Otherwise we can distinguish between two cases. First \eqn{Lambda} is not specified, hence the list is four-dimensional. The first element determines the CutOff value of the data by using the sorted CutOff vector (using the plot, this means the row in which the graph is shown). The second element specifies the maximal number of modes by using the sorted mode vector (again using the plot, this means the column of the plot). The third element selects the \eqn{\lambda} of the graph. For each plot and each \eqn{\lambda}, the following information is stored: the value of \eqn{\lambda}, the \eqn{\lambda}-clusters and the excess mass vector. Using the default setting \eqn{[2,2,5,2]} shows the \eqn{\lambda}-clusters of the fifth smallest \eqn{\lambda} of the \eqn{CutOff=2-M=2}-plot.

If \eqn{Lambda} is declared manually, the list is three-dimensional. Hence, the first argument denotes the maximal number of modes (the column of the graph). The second argument indicates the \eqn{\lambda}  by the position held by it in the \eqn{Lambda} vector. As in "exmplot" only two information are shown. The \eqn{\lambda}-clusters (\eqn{[,,1]}) and the vector of excess mass (\eqn{[,,2]}), as the value of \eqn{\lambda} is known. 
}
\references{Muller, D. W. and Sawitzki, G., 09.1991, "Excess Mass Estimates and Tests for Multimodality", \emph{Journal of the American Statistical Association} \bold{, Vol. 86, No. 415, pp. 738--746}, \url{http://www.jstor.org/stable/2290406}
}
\author{Marc-Daniel Mildenberger \email{mildenberger.stat@web.de}, based on earlier code from Dr. Guenther Sawitzki \email{gs@statlab.uni-heidelberg.de} 
}

\seealso{\code{\link{excessm}, \link{exmplot}, \link{exmsilhouette}}
}
\examples{
library(MASS)
attach(geyser)

##calculating excess mass plots for duration of 'Old Faithful Geyser',
##specifying CutOff and number of steps manually
mexmsilhouette(duration, CutOff=c(1,2), steps=60)

##Allowing for three different maximal number of modes 
##and CutOff factors as well as color. 
##The rug plot is omitted and numerical data is requested.
res=mexmsilhouette(duration, M=c(2,3,7), CutOff=c(0.8,1,2), col=TRUE, rug=FALSE, rdata=TRUE)

##Lambda is specified, color is set to true, numerical data is requested
L=seq(.01,.25,0.005)
res=mexmsilhouette(duration, M=c(2,3,4), Lambda=L, col=TRUE, rdata=TRUE)
}

\keyword{ hplot }
\keyword{ htest }
\keyword{ univ }
\keyword{ nonparametric }
