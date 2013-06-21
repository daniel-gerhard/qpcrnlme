\name{print-methods}
\alias{print.ddct}
\alias{print.nlmect}
\alias{coef.nlmect}
\alias{coefficients.nlmect}
\alias{vcov.nlmect}

\title{nlmect and ddct methods}
\description{print and other methods for nlmect and ddct objects.}
\usage{
\method{print}{ddct}(x, digits = 3, ...)
\method{print}{nlmect}(x, digits = 3, ...)
\method{coef}{nlmect}(object, ...)
\method{coefficients}{nlmect}(object, ...)
\method{vcov}{nlmect}(object, ...)
}
\arguments{
  \item{x}{An object of class nlmect or ddct}
  \item{object}{An object of class nlmect}
  \item{digits}{number of digits}
  \item{\dots}{\dots}
}

\author{Daniel Gerhard}

\seealso{\code{\link{qpcr_nlme}}, \code{\link{qpcr_nlme_formula}}, \code{\link{ddctcomp}}}

\keyword{methods}
