\name{qpcr_nlme}
\alias{qpcr_nlme}
\title{5-parameter logistic mixed model analysis of qpcr data}
\description{Fits a nonlinear mixed model to qpcr fluorescence data,
  allowing for one factorial ANOVA designs with several levels for a
  single treatment factor at multiple genes. Marginal c(t) values are calculated.}
\usage{
qpcr_nlme(response, cycle, gene, treatment, brep, well, data, 
 cutoff, nGQ=5, verbose=TRUE)
}
\arguments{
  \item{response}{Name of a column in data, denoting a numeric response
    vector with fluorescence observations.}
  \item{cycle}{Name of a column in data, denoting a numeric response
    vector with cycle numbers.}
  \item{gene}{Name of a column in data, denoting a factor coding for
    different genes.}
  \item{treatment}{Name of a column in data, denoting a factor coding for
    several treatment groups}
  \item{brep}{Name of a column in data, denoting a factor coding for
    biological replications.}
  \item{well}{Name of a column in data, denoting a factor coding for
    technical replications.}
  \item{data}{A data.frame object.}
  \item{cutoff}{cutoff value t to define c(t) cycle thresholds.}
  \item{nGQ}{Number of nodes and weights for Gaussian Quadrature}
  \item{verbose}{Some text output during calculation...}
}
\value{An object of class ddct.}

\author{Daniel Gerhard <gerhard@biostat.uni-hannover.de>}


\seealso{\code{\link{nlme}}}

\examples{
\dontrun{

data(pt6c)
pt6c$Time <- as.factor(pt6c$Time)
pt6result <- qpcr_nlme(response="Fluorescence", cycle="Cycle", gene="Target", 
                       treatment="Time", brep="Content", well="Well", 
                       data=pt6c, cutoff=100, verbose=TRUE)

data(pt3c)
pt3c$Time <- as.factor(pt3c$Time)
pt3result <- qpcr_nlme(response="Fluorescence", cycle="Cycle", gene="Target", 
                       treatment="Time", brep="Content", well="Well", 
                       data=pt3c, cutoff=100, verbose=TRUE)

data(KOTWTcurves)
Kresult <- qpcr_nlme(response="Fluorescence", cycle="Cycle", gene="Target", treatment="Sample", brep="Content", well="Well", 
                     data=KOTWTcurves, cutoff=100, verbose=TRUE)

data(QTTWTcurves)
Qresult <- qpcr_nlme(response="Fluorescence", cycle="Cycle", gene="Target", treatment="Sample", brep="Content", well="Well", 
                     data=QTTWTcurves, cutoff=10, nGQ=5, verbose=TRUE)
}
}


\keyword{ misc }
