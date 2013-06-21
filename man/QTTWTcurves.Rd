\name{QTTWTcurves}
\alias{QTTWTcurves}
\docType{data}
\title{QTTWTcurves}
\description{Dataset from the Institute of plant nutrition (LUH). In the
  experiments an insertion has introduced into the genome of different
  rice lines. This insertion changes gene expression from the gene of
  interest. By using the quantitative RT-PCR the changes in gene
  expression of a target gene will analyse between different samples and
  the wild type. 
  Content: ID of biological replicates
  Target: gene of interest and EF (elongation factor) gene of control
  Samples: rice lines with an insertion (different letters) and wild type (W)}
\usage{data(QTTWTcurves)}
\format{
  A data frame with 3780 observations on the following 6 variables.
  \describe{
    \item{\code{Well}}{a factor with well IDs}
    \item{\code{Content}}{factor levels denoting biological replicates}
    \item{\code{Target}}{2 different genes 10g30610 and EF (control)}
    \item{\code{Sample}}{treatments Q, TT, WT}   
    \item{\code{Cycle}}{cycle number} 
    \item{\code{Fluorescence}}{observed fluorescence (proportional to PCR material)}
  }
}
\author{Eline Biedermann <eline.biedermann@gmx.de>}

\examples{
data(QTTWTcurves)
ggplot(QTTWTcurves, aes(x=Cycle, y=Fluorescence, colour=Sample,
group=Well)) + geom_line() + facet_wrap(~Target)

\dontrun{
qttwt <- qpcr_nlme(response="Fluorescence", cycle="Cycle", gene="Target", 
                   treatment="Sample", brep="Content", well="Well", 
                   data=QTTWTcurves, cutoff=10, nGQ=5, verbose=TRUE)
ddctcomp(qttwt, control_gene="EF", control_treatment=NULL, 
         conf.level=0.95, adjusted=TRUE, ratio_ddct=TRUE)
}

}

\keyword{datasets}
