\name{ExcessMass}
\alias{ExcessMass-package}
\alias{ExcessMass}
\docType{package}
\title{Excess Mass Calculation and Plots
}
\description{
Implementation of a function which calculates the empirical excess mass for a given \eqn{\lambda} and given maximal number of modes (excessm). Offering powerful plot features to visualize empirical excess mass (exmsilhouette). This includes the possibility of drawing several plots (with different maximal number of modes / cut off values) in a single graph. Furthermore, plotting the empirical excess mass against lambda is implemented (exmplot).
}
\details{
\tabular{ll}{
Package: \tab ExcessMass\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2017-05-17\cr
License: \tab GPL
}
}
\author{Marc-Daniel Mildenberger \email{mildenberger.stat@web.de}, based on earlier code from Dr. Guenther Sawitzki \email{gs@statlab.uni-heidelberg.de}
}
\references{Muller, D. W. and Sawitzki, G., 09.1991, "Excess Mass Estimates and Tests for Multimodality", \emph{Journal of the American Statistical Association} \bold{, Vol. 86, No. 415, pp. 738--746}, \url{http://www.jstor.org/stable/2290406}

Muller, D. W., 12.1992, "The Excess Mass Approach in Statistics", \emph{Beitraege zur Statistik -- StatLab Heidelberg}, \url{http://archiv.ub.uni-heidelberg.de/volltextserver/21357/}}
\keyword{ package }

\examples{
library(MASS)
attach(geyser)
excessm(duration, lambda=0.2)

x <- rnorm(1000)
exmsilhouette(x, M=2, CutOff=0.5)

mexmsilhouette(duration, CutOff=c(1,2), steps=60)
}
