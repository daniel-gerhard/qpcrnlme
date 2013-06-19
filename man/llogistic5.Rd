\name{llogistic5}
\alias{llogistic5}
\title{5-parameter log-logistic model function}
\description{The general asymmetric five-parameter log-logistic model for describing dose-response relationships.}
\usage{
llogistic5(cycle, b, c, d, e, f)
}
\arguments{
  \item{cycle}{cycle number}
  \item{b}{slope}
  \item{c}{lower asymptote}
  \item{d}{upper asymptote}
  \item{e}{cycle number at inflection point}
  \item{f}{asymmetry parameter}
}

\author{Daniel Gerhard <gerhard@biostat.uni-hannover.de>}

\seealso{\code{\link{nlme}}}

\keyword{ misc }
