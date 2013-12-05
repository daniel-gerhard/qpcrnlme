ddctcomp <- function(object, control_gene=NULL, control_treatment=NULL, conf.level=0.95, adjusted=FALSE, ratio_ddct=TRUE){
  nt <- object$nt
  ng <- object$ng
  est <- object$ct[,1]
  vce <- object$vcov
  df <- object$df
  dat <- object$dat
  if (ratio_ddct == TRUE) dsymb <- " / " else dsymb <- " - "
  if (!is.null(control_gene)) if (!(control_gene %in% levels(dat$gene))) stop(paste(control_gene, "is not a gene name!"))
  if (!is.null(control_treatment)) if (!(control_treatment %in% levels(dat$treatment))) stop(paste(control_treatment, "is not a treatment level!"))
  ##### ddct comparisons
  # define pairwise comparisons
  # between treatments
  if (is.null(control_treatment)){
    tcombi <- combn(1:nt, 2)
  } else {
    wtr <- which(levels(dat$treatment) == control_treatment)
    tcombi <- rbind(wtr, (1:nt)[-wtr])
  }
  # between genes
  if (is.null(control_gene)){
    gcombi <- combn(0:(ng-1), 2)
  } else {
    wge <- which(levels(dat$gene) == control_gene)
    gcombi <- rbind(wge, (1:ng)[-wge])-1
  }
  
  # estimate
  ddctest <- as.vector(apply(gcombi, 2, function(gc){
    apply(tcombi, 2, function(tc){
      if (ratio_ddct == TRUE){
        (est[tc[1] + nt*gc[1]] /  est[tc[2] + nt*gc[1]]) / (est[tc[1] + nt*gc[2]] / est[tc[2] + nt*gc[2]])
      } else {
        (est[tc[1] + nt*gc[1]] -  est[tc[2] + nt*gc[1]]) - (est[tc[1] + nt*gc[2]] - est[tc[2] + nt*gc[2]])
      }
    })
  }))
  
  # variance estimate
  forml <- apply(gcombi, 2, function(gc){
    apply(tcombi, 2, function(tc){
      as.formula(paste("~ ", "(", "x", tc[1] + nt*gc[1], dsymb, "x", tc[2] + nt*gc[1], ")", dsymb, "(", "x", tc[1] + nt*gc[2], dsymb, "x", tc[2] + nt*gc[2], ")", sep=""))
    })
  })[[1]]
  syms <- paste("x", 1:length(est), sep = "")
  for (i in 1:length(est)) assign(syms[i], est[i])
  gec <- t(sapply(forml, function(form) as.numeric(attr(eval(deriv(form, syms)), "gradient"))))
  vcc <- gec %*% vce %*% t(gec)
  ddctstd <- sqrt(diag(vcc))
  
  # ddct names
  lt <- abbreviate(levels(dat$treatment),3)
  lg <- abbreviate(levels(dat$gene),3)
  nm <- apply(gcombi+1, 2, function(gc){
    apply(tcombi, 2, function(tc){
      paste(paste("(",paste(paste(lt[tc[1]], lg[gc[1]], sep=":"), paste(lt[tc[2]], lg[gc[1]], sep=":"), sep=dsymb), ")", sep=""), paste("(",paste(paste(lt[tc[1]], lg[gc[2]], sep=":"), paste(lt[tc[2]], lg[gc[2]], sep=":"), sep=dsymb), ")", sep=""), sep=dsymb)
    })    
  })
  
  
  # confidence intervals
  if (adjusted == TRUE){
    cr <- cov2cor(vcc)
    quant <- qmvt(conf.level, tail="both.tails", df=df, corr=cr)$quantile
  } else {
    quant <- qt(0.5+0.5*conf.level, df=df)
  }
  
  ## p-values
  if (ratio_ddct == TRUE) deltashift <- 1 else deltashift <- 0
  tstat <- (ddctest-deltashift)/ddctstd
  if (adjusted == TRUE){
    cr <- cov2cor(vcc)
    dim <- ncol(cr)
    pfct <- function(q){
      low <- rep(-abs(q), dim)
      upp <- rep(abs(q), dim)
      pmvt(lower = low, upper = upp, df = df, corr = cr)
    }
    pv <- numeric(length(tstat))
    for (i in 1:length(tstat)) pv[i] <- 1 - pfct(tstat[i])
  } else {
    pv <- pmin(1, 2*pt(abs(tstat), df=df, lower.tail=FALSE)*2)
  }
  
  ci <- cbind(ddctest, ddctstd, ddctest-quant*ddctstd, ddctest+quant*ddctstd, pv)
  rownames(ci) <- nm
  colnames(ci) <- c("estimate", "std.err", "lower", "upper","p-value")
  
  out <- list()
  out$vcovddct <- vcc
  out$coefmat <- ci
  out$df <- df
  class(out) <- "ddct"
  return(out)
}
