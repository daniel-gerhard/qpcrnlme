print.nlmect <- function(x, digits = 3, ...){
  cat("\nc(t) estimates:\n")
  print(x$ct, digits=digits)
}

print.ddct <- function(x, digits = 3, ...){
  cat("\ndelta delta c(t) estimates:\n")
  printCoefmat(x$coefmat, digits=digits, P.values=TRUE, has.Pvalue=TRUE)
}

coef.nlmect <- function(object, ...){
  return(object$ct[,1])
}
coefficients.nlmect <- function(object, ...){
  return(object$ct[,1])
}

vcov.nlmect <- function(object, ...){
  vc <- object$vcov
  rownames(vc) <- colnames(vc) <- names(object$ct[,1])
  return(vc)
}
