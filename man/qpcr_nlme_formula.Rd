\name{qpcr_nlme_formula}
\alias{qpcr_nlme_formula}
\title{5-parameter logistic mixed model analysis of qpcr data}
\description{Fits a nonlinear mixed model to qpcr fluorescence data,
  allowing for one factorial ANOVA designs with several levels for a
  single treatment factor at multiple genes. Marginal c(t) values and
  corresponding delta delta c(t) values are calculated.}
\usage{
qpcr_nlme_formula(response, cycle, gene, trtformula, brep, well, data, newdata, cutoff, nGQ=5, verbose=TRUE)
}
\arguments{
  \item{response}{Name of a column in data, denoting a numeric response
    vector with fluorescence observations.}
  \item{cycle}{Name of a column in data, denoting a numeric response
    vector with cycle numbers.}
  \item{gene}{Name of a column in data, denoting a factor coding for
    different genes.}
  \item{trtformula}{One-sided formula to model treatment effects}
  \item{brep}{Name of a column in data, denoting a factor coding for
    biological replications.}
  \item{well}{Name of a column in data, denoting a factor coding for
    technical replications.}
  \item{data}{A data.frame object.}
  \item{newdata}{A data.frame object with covariate values at which
    cycle thresholds are estimated.}
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
qpcr_nlme_formula(response="Fluorescence", cycle="Cycle", gene="Target", trtformula=~ poly(Time, 4, raw=TRUE), brep="Content", well="Well", data=pt6c, newdata=data.frame(Time=seq(1,8, length=5)), cutoff=100, nGQ=5, verbose=TRUE)

}
}


\keyword{ misc }
