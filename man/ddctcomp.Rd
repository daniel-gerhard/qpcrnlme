\name{ddctcomp}
\alias{ddctcomp}
\title{Delta Delta c(t) calculation}
\description{All-pair comparisons of treatment levels or comparisons to a control treatment by Delta Delta c(t) parameters given a nlmect object.}
\usage{
ddctcomp(object, control_gene=NULL, control_treatment=NULL, 
         conf.level=0.95, adjusted=FALSE, ratio_ddct=TRUE)
}
\arguments{
  \item{object}{An object of class nlmect obtained by calling the function qpcr_nlme.}
  \item{control_gene}{Name of a single control gene in the vector of
    gene names.}
  \item{control_treatment}{Name of a single control treatment group in
    the vector of treatment names.}
  \item{conf.level}{confidence level for delta delta c(t) confidence
    intervals}
  \item{adjusted}{if TRUE, simultaneous confidence intervals controlling
  the FWER are reported.}
  \item{ratio_ddct}{If TRUE delta delta c(t) is defined as ratio,
    otherwise the difference of c(t) values}
}
\value{An object of class ddct.}

\author{Daniel Gerhard <gerhard@biostat.uni-hannover.de>}


\seealso{\code{\link{nlme}}}

\examples{
\dontrun{

data(pt6c)
pt6c$Time <- as.factor(pt6c$Time)
pt6result <- qpcr_nlme(response="Fluorescence", cycle="Cycle", gene="Target", treatment="Time", 
                       brep="Content", well="Well", data=pt6c, cutoff=100, verbose=TRUE)
ddctcomp(pt6result, control_gene="eEf", control_treatment="0", 
         conf.level=0.95, adjusted=TRUE, ratio_ddct=TRUE)

}
}


\keyword{ misc }
