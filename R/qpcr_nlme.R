qpcr_nlme <- function(response, cycle, gene, control_gene=NULL, treatment, control_treatment=NULL, brep, well, data, cutoff, conf.level=0.95, adjusted=FALSE, ratio_ddct=TRUE, nGQ=5, verbose=TRUE){
  require(drc)
  require(nlme)
  require(statmod)

  # function argument checks
  if (!("data.frame" %in% class(data))) stop("data is not of class data.frame")
  dnames <- names(data)
  if (!(response %in% dnames)) stop(paste("response=", response, "is not a variable in the dataframe!"))
  if (!(cycle %in% dnames)) stop(paste("cycle=", cycle, "is not a variable in the dataframe!"))
  if (!(gene %in% dnames)) stop(paste("gene=", gene, "is not a variable in the dataframe!"))
  if (!(treatment %in% dnames)) stop(paste("treatment=", treatment, "is not a variable in the dataframe!"))
  if (!(brep %in% dnames)) stop(paste("brep=", brep, "is not a variable in the dataframe!"))
  if (!(well %in% dnames)) stop(paste("well=", well, "is not a variable in the dataframe!"))
  if (!is.factor(data[,gene])) stop(paste(gene, "is not a factor!"))
  if (!is.factor(data[,treatment])) stop(paste(treatment, "is not a factor!"))
  if (length(levels(data[,gene])) < 2) stop(paste(gene, "has less than 2 levels!"))
  if (length(levels(data[,treatment])) < 2) stop(paste(treatment, "has less than 2 levels!"))
  if (!is.null(control_gene)) if (!(control_gene %in% levels(data[,gene]))) stop(paste(control_gene, "is not a level in", gene, "!"))
  if (!is.null(control_treatment)) if (!(control_treatment %in% levels(data[,treatment]))) stop(paste(control_treatment, "is not a level in", treatment, "!"))

  # working dataset
  dat <- data.frame(response=data[,response], cycle=data[,cycle],gene=data[,gene],treatment=data[,treatment],brep=data[,brep],well=data[,well])
  ng <- length(levels(dat$gene))
  nt <- length(levels(dat$treatment))
  dat$gt <- as.factor(paste(dat$gene, dat$treatment, sep=":"))
  if (ratio_ddct == TRUE) dsymb <- " / " else dsymb <- " - "
  
  # starting values
  if (verbose) cat("Searching for starting values ... ")
  dfm <- drm(response ~ cycle, curveid=gt, data=dat, fct=LL.5(), separate=TRUE)
  cf <- matrix(coefficients(dfm), ncol=5, byrow=TRUE)
  start <- c(as.vector(cf[,c(1,2,4,5)]), mean(cf[,3]))
  if (verbose) cat("done\n")
  
  # nonlinear mixed model
  if (verbose) cat("Parameter estimation ... ")
  fm <- nlme(response ~ llogistic5(cycle, b, c, d, e, f),
           fixed = list(b + c + e + f ~ gt-1, d ~ 1),
           random = list(brep=pdDiag(d + e + b ~ gene-1), well=pdDiag(d + e + b ~ 1)), 
           start=start, data=dat, control=nlmeControl(maxIter = 500), method="ML")

  if (verbose) cat("done\n")
  fix <- fixef(fm)
  vc <- VarCorr(fm)
  vcf <- vcov(fm)

  ### marginal ct calculation  
  if (verbose) cat("Calculating marginal c(t) ... ")
  rsd <- suppressWarnings(as.numeric(vc))

  # extract variance-components
  oew <- length(rsd)-2
  obw <- length(rsd)-1
  oe <- ((length(rsd)/2) + 1 + ng + 1:ng)
  ob <- ((length(rsd)/2) + 1 + 2*ng + 1:ng)

  # compute Gaussian-Quadrature nodes and weights (package statmod)
  oute <- lapply(oe, function(i) gauss.quad.prob(nGQ,dist="normal",mu=0,sigma=rsd[i]))  
  outb <- lapply(ob, function(i) gauss.quad.prob(nGQ,dist="normal",mu=0,sigma=rsd[i]))
  outew <- gauss.quad.prob(nGQ,dist="normal",mu=0,sigma=rsd[oew])
  outbw <- gauss.quad.prob(nGQ,dist="normal",mu=0,sigma=rsd[obw])


  # extract fixed effects
  fixm <- matrix(fix[-length(fix)], ncol=4)
  fb <- fixm[,1]
  fc <- fixm[,2]
  fd <- fix[length(fix)]
  fe <- fixm[,3]
  ff <- fixm[,4]

  if (cutoff < max(fc) | cutoff > fd) stop(paste("cutoff=", cutoff, "is outside the limits [", round(max(fc),2), ";", round(fd,2), "]."))

  # marginal c(t) function
  invreg <- function(cutoff, fb, fc, fe, ff, fd, outb, oute, outbw, outew, ng, nt){
    yopt <- function(x, fx, cutoff, outb, oute, outbw, outew){
      sapply(x, function(x) (sum(outb$weights*sapply(outb$nodes, function(bn)
          sum(outbw$weights*sapply(outbw$nodes, function(bw)
          sum(oute$weights*sapply(oute$nodes, function(en)
          sum(outew$weights*logistic5(x,
                                      fx[1] + bw + bn,
                                      fx[2],
                                      fx[5],
                                      fx[3] + en + outew$nodes,
                                      fx[4])))))))) - cutoff))
    }
    fbm <- matrix(fb, ncol=ng)
    fcm <- matrix(fc, ncol=ng)
    fem <- matrix(fe, ncol=ng)
    ffm <- matrix(ff, ncol=ng)
    # Ct estimation
    cts <- sapply(1:ng, function(j){
      sapply(1:nt, function(i){
        uniroot(yopt, interval=c(-cutoff,fd),
                fx=c(fbm[i,j], fcm[i,j], fem[i,j], ffm[i,j], fd),
                cutoff=cutoff,
                outb=outb[[j]], oute=oute[[j]], outbw=outbw, outew=outew)$root
      })
    })
    return(as.vector(cts))
  }

  # c(t) estimation and gradient calculation
  est <- numericDeriv(quote(invreg(cutoff, fb, fc, fe, ff, fd, outb, oute, outbw, outew, ng, nt)), c("fb","fc","fe","ff","fd"))
  # Delta-method for variance estimation
  vce <- attr(est,"gradient") %*% vcf %*% t(attr(est,"gradient"))
  std <- sqrt(diag(vce))

  ctest <- cbind("estimate"=est, "std.err"=std)
  rownames(ctest) <- levels(dat$gt)
  if (verbose) cat("done\n")
  
  ##### ddct comparisons
  # define pairwise comparisons
  # between treatments
  if (verbose) cat("Calculating delta delta c(t) ... ")
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
  if (verbose) cat("done\n")
  
  ## df  (using lme containment df for aggregated data)
  if (verbose) cat("Computing degrees of freedom  ... ")
  fms <-  drm(response ~ cycle, curveid=well, data=dat, fct=L.5(), separate=TRUE)
  fcs <- matrix(coefficients(fms), ncol=5, byrow=TRUE)
  sct <- apply(fcs, 1, function(fc) ((log( ((fc[3]-fc[2])/( cutoff-fc[2]))^(1/fc[5]) - 1) / fc[1]) + fc[4]  ))
  spdat <- data.frame(treatment=levels(dat$treatment)[tapply(dat$treatment, dat$well, unique)],
                      gene=levels(dat$gene)[tapply(dat$gene, dat$well, unique)],
                      ct=sct, brep=tapply(dat$brep, dat$well, unique))
  aggdat <- aggregate(ct ~ treatment + gene + brep, data=na.omit(spdat), mean)
  fmldf <- lme(ct ~ treatment*gene, na.omit(aggdat), random=~gene-1|brep)
  df <- anova(fmldf)$denDF[4]
  if (verbose) cat("done\n")

  # ddct names
  lt <- abbreviate(levels(dat$treatment),3)
  lg <- abbreviate(levels(dat$gene),3)
  nm <- apply(gcombi+1, 2, function(gc){
    apply(tcombi, 2, function(tc){
      paste(paste("(",paste(paste(lt[tc[1]], lg[gc[1]], sep=":"), paste(lt[tc[2]], lg[gc[1]], sep=":"), sep=dsymb), ")", sep=""), paste("(",paste(paste(lt[tc[1]], lg[gc[2]], sep=":"), paste(lt[tc[2]], lg[gc[2]], sep=":"), sep=dsymb), ")", sep=""), sep=dsymb)
    })    
  })
  

  # confidence intervals
  if (verbose) cat("Estimating confidence limits  ... ")
  if (adjusted == TRUE){
    require(mvtnorm)
    cr <- cov2cor(vcc)
    quant <- qmvt(conf.level, tail="both.tails", df=df, corr=cr)$quantile
  } else {
    quant <- qt(0.5+0.5*conf.level, df=df)
  }

  ## p-values
  if (ratio_ddct == TRUE) deltashift <- 1 else deltashift <- 0
  tstat <- (ddctest-deltashift)/ddctstd
  if (adjusted == TRUE){
    require(mvtnorm)
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
  if (verbose) cat("done\n")
  
  out <- list()
  out$nlme <- fm
  out$ct <- ctest
  out$vcov <- vce
  out$vcovddct <- vcc
  out$coefmat <- ci
  out$df <- df
  out$gradient <- attr(est, "gradient")
  class(out) <- "ddct"
  return(out)
}


print.ddct <- function(x, digits = 3, ...){
  cat("\nc(t) estimates:\n")
  print(x$ct, digits=digits)
  cat("\ndelta delta c(t) estimates:\n")
  printCoefmat(x$coefmat, digits=digits, P.values=TRUE, has.Pvalue=TRUE)
}

