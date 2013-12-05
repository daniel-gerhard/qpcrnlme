llogistic5 <- function(cycle, b, c, d, e, f){
  .value <- c + (d - c)/((1 + exp(b * (log(cycle) - log(e))))^f)
  .actualArgs <- as.list(match.call()[c("b", "c", "d", "e", "f")])
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    dfct <- function (cycle, parm){
      parmMat <- matrix(c(NA, NA, NA, NA, NA), nrow(parm), 5, byrow = TRUE)
      parmMat[,] <- parm
      t1 <- parmMat[, 3] - parmMat[, 2]
      t2 <- exp(parmMat[, 1] * (log(cycle) - log(parmMat[, 4])))
      t5 <- (1 + t2)^parmMat[, 5]
      eval(parse(text="cbind(-t1 * drc:::xlogx(cycle/parmMat[, 4], parmMat[, 1], parmMat[,5] + 1) * parmMat[, 5], 1 - 1/t5, 1/t5, t1 * parmMat[,5] * drc:::divAtInf(t2, (1 + t2)^(parmMat[, 5] + 1)) * parmMat[,1]/parmMat[, 4], -t1 * drc:::divAtInf(log(1 + t2), t5))"))     
    }    
    .grad <- dfct(cycle, cbind(b, c, d, e, f))
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}


logistic5 <- function(cycle, b, c, d, e, f){
  .value <- c + (d - c)/((1 + exp(b * (cycle - e)))^f)
  .actualArgs <- as.list(match.call()[c("b", "c", "d", "e", "f")])
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    dfct <- function (cycle, parm){
      parmMat <- matrix(c(NA, NA, NA, NA, NA), nrow(parm), 5, byrow = TRUE)
      parmMat[,] <- parm
      t1 <- parmMat[, 3] - parmMat[, 2]
      t2 <- exp(parmMat[, 1] * (cycle - parmMat[, 4]))
      t3 <- (1 + t2)^(2 * parmMat[, 5])
      t4 <- parmMat[, 5] * ((1 + t2)^(-parmMat[, 5] - 1))
      t5 <- (1 + t2)^(parmMat[, 5])
      cbind(-t1 * t2 * t4 * (cycle - parmMat[, 4]), 1 - 1/t5, 1/t5, t1 * t2 * t4 * parmMat[, 1], -t1 * log(1 + t2)/t5)
    }    
    .grad <- dfct(cycle, cbind(b, c, d, e, f))
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

