\name{pt3c}
\alias{pt3c}
\docType{data}
\title{Fluorescence curves for PT3 expression at phosphate deficiency}
\description{RT-PCR analysis of the phosphate transporter PT3 in rice at
  several time periods of phosphate deficiency.}
\usage{data(pt3c)}
\format{
  A data frame with 3600 observations on the following 6 variables.
  \describe{
    \item{\code{Cycle}}{cycle number}
    \item{\code{Fluorescence}}{observed fluorescence (proportional to PCR material)}
    \item{\code{Well}}{a factor with well IDs}
    \item{\code{Content}}{factor levels denoting biological replicates}
    \item{\code{Time}}{0,2,4,6, and 8 hours in a nutrient solution at
      low phosphate content}
    \item{\code{Target}}{2 different genes PT3 and eEf (control)}
  }
}
\author{Melanie Bremer <bremer@pflern.uni-hannover.de>}


\examples{
data(pt3c)
ggplot(pt3c, aes(x=Cycle, y=Fluorescence, colour=Time, group=Well)) +
geom_line()

\dontrun{
pt3c$Time <- as.factor(pt3c$Time)
pt3result <- qpcr_nlme(response="Fluorescence", cycle="Cycle",
gene="Target", control_gene="eEf", treatment="Time",
control_treatment="0", brep="Content", well="Well", data=pt3c,
cutoff=100, conf.level=0.95, adjusted=TRUE, ratio_ddct=TRUE,
verbose=TRUE)
}


}


\keyword{datasets}
