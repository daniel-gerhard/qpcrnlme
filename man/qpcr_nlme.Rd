\name{qpcr_nlme}
\alias{qpcr_nlme}
\title{5-parameter logistic mixed model analysis of qpcr data}
\description{Fits a nonlinear mixed model to qpcr fluorescence data,
  allowing for one factorial ANOVA designs with several levels for a
  single treatment factor at multiple genes. Marginal c(t) values and
  corresponding delta delta c(t) values are calculated.}
\usage{
qpcr_nlme(response, cycle, gene, control_gene=NULL, treatment,
control_treatment=NULL, brep, well, data, cutoff, conf.level=0.95,
adjusted=FALSE, ratio_ddct=TRUE, nGQ=5, verbose=TRUE)
}
\arguments{
  \item{response}{Name of a column in data, denoting a numeric response
    vector with fluorescence observations.}
  \item{cycle}{Name of a column in data, denoting a numeric response
    vector with cycle numbers.}
  \item{gene}{Name of a column in data, denoting a factor coding for
    different genes.}
  \item{control_gene}{Name of a single control gene in the vector of
    gene names.}
  \item{treatment}{Name of a column in data, denoting a factor coding for
    several treatment groups}
  \item{control_treatment}{Name of a single control treatment group in
    the vector of treatment names.}
  \item{brep}{Name of a column in data, denoting a factor coding for
    biological replications.}
  \item{well}{Name of a column in data, denoting a factor coding for
    technical replications.}
  \item{data}{A data.frame object.}
  \item{cutoff}{cutoff value t to define c(t) cycle thresholds.}
  \item{conf.level}{confidence level for delta delta c(t) confidence
    intervals}
  \item{adjusted}{if TRUE, simultaneous confidence intervals controlling
  the FWER are reported.}
  \item{ratio_ddct}{If TRUE delta delta c(t) is defined as ratio,
    otherwise the difference of c(t) values}
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
pt6result <- qpcr_nlme(response="Fluorescence", cycle="Cycle", gene="Target", control_gene="eEf", treatment="Time", control_treatment="0", brep="Content", well="Well", data=pt6c, cutoff=100, conf.level=0.95, adjusted=TRUE, ratio_ddct=TRUE, verbose=TRUE)

data(pt3c)
pt3c$Time <- as.factor(pt3c$Time)
pt3result <- qpcr_nlme(response="Fluorescence", cycle="Cycle", gene="Target", control_gene="eEf", treatment="Time", control_treatment="0", brep="Content", well="Well", data=pt3c, cutoff=100, conf.level=0.95, adjusted=TRUE, ratio_ddct=TRUE, verbose=TRUE)

data(KOTWTcurves)
Kresult <- qpcr_nlme(response="Fluorescence", cycle="Cycle", gene="Target", control_gene="EF", treatment="Sample", control_treatment=NULL, brep="Content", well="Well", data=KOTWTcurves, cutoff=100, conf.level=0.95, adjusted=TRUE, ratio_ddct=TRUE, verbose=TRUE)

data(QTTWTcurves)
Qresult <- qpcr_nlme(response="Fluorescence", cycle="Cycle", gene="Target", control_gene="EF", treatment="Sample", control_treatment=NULL, brep="Content", well="Well", data=QTTWTcurves, cutoff=10, conf.level=0.95, adjusted=TRUE, ratio_ddct=TRUE, nGQ=5, verbose=TRUE)
}
}


\keyword{ misc }
