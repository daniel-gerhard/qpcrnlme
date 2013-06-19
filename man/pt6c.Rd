\name{pt6c}
\alias{pt6c}
\docType{data}
\title{Fluorescence curves for PT6 expression at phosphate deficiency}
\description{RT-PCR analysis of the phosphate transporter PT6 in rice at
  several time periods of phosphate deficiency.}
\usage{data(pt6c)}
\format{
  A data frame with 3600 observations on the following 6 variables.
  \describe{
    \item{\code{Cycle}}{cycle number}
    \item{\code{Fluorescence}}{observed fluorescence (proportional to PCR material)}
    \item{\code{Well}}{a factor with well IDs}
    \item{\code{Content}}{factor levels denoting biological replicates}
    \item{\code{Time}}{0,2,4,6, and 8 hours in a nutrient solution at
      low phosphate content}
    \item{\code{Target}}{2 different genes PT6 and eEf (control)}
  }
}
\author{Melanie Bremer <bremer@pflern.uni-hannover.de>}

\examples{
data(pt6c)
ggplot(pt6c, aes(x=Cycle, y=Fluorescence, colour=Time, group=Well)) +
geom_line()

\dontrun{

# definition of delta delta c(t)
library(msm) # contains deltamethod function
p <- 25 # number of time points for prediction
# (computation time increases with increasing p)
tt <- seq(0,8, length=p) # time range for prediction 
forms <- sapply(2:p, function(i){
   as.formula(paste("~ (x",i," / x",i+p,") / (x",1," / x",p+1,")", sep="")))
}

## estimation of marginal cycle thresholds c(t)
# 1 degree polynomial
cte1 <- qpcr_nlme_formula(response="Fluorescence", cycle="Cycle", gene="Target",
                          trtformula=~poly(Time, 1, raw=TRUE),
                          brep="Content", well="Well",
                          data=pt6c, newdata=data.frame(Time=tt),
                          cutoff=100, nGQ=5, verbose=TRUE)
sds1 <- deltamethod(forms, cte1$ct[,1], cte1$vcov)
est1 <- sapply(2:p, function(i) (cte1$ct[i,1] / cte1$ct[i+p,1]) / (cte1$ct[1,1] / cte1$ct[p+1,1]))
# 2 degree polynomial
cte2 <- qpcr_nlme_formula(response="Fluorescence", cycle="Cycle", gene="Target",
                          trtformula=~poly(Time, 2, raw=TRUE),
                          brep="Content", well="Well",
                          data=pt6c, newdata=data.frame(Time=tt),
                          cutoff=100, nGQ=5, verbose=TRUE)
sds2 <- deltamethod(forms, cte2$ct[,1], cte2$vcov)
est2 <- sapply(2:p, function(i) (cte2$ct[i,1] / cte2$ct[i+p,1]) / (cte2$ct[1,1] / cte2$ct[p+1,1]))
# 3 degree polynomial
cte3 <- qpcr_nlme_formula(response="Fluorescence", cycle="Cycle", gene="Target",
                          trtformula=~poly(Time, 3, raw=TRUE),
                          brep="Content", well="Well",
                          data=pt6c, newdata=data.frame(Time=tt),
                          cutoff=100, nGQ=5, verbose=TRUE)
sds3 <- deltamethod(forms, cte3$ct[,1], cte3$vcov)
est3 <- sapply(2:p, function(i) (cte3$ct[i,1] / cte3$ct[i+p,1]) / (cte3$ct[1,1] / cte3$ct[p+1,1]))
# 4 degree polynomial
cte4 <- qpcr_nlme_formula(response="Fluorescence", cycle="Cycle", gene="Target",
                          trtformula=~poly(Time, 4, raw=TRUE),
                          brep="Content", well="Well",
                          data=pt6c, newdata=data.frame(Time=tt),
                          cutoff=100, nGQ=5, verbose=TRUE)
sds4 <- deltamethod(forms, cte4$ct[,1], cte4$vcov)
est4 <- sapply(2:p, function(i) (cte4$ct[i,1] / cte4$ct[i+p,1]) / (cte4$ct[1,1] / cte4$ct[p+1,1]))


## model weights based on information criteria
ic <- BIC(cte1$nlme, cte2$nlme, cte3$nlme, cte4$nlme) 
dm <- ic$BIC - min(ic$BIC)
w <- exp(-0.5*dm)/sum(exp(-0.5*dm))
round(w,3)


### summary of estimation results
edat <- data.frame(time=rep(tt[-1], 4),
                   est=c(est1,est2,est3,est4),
                   std=c(sds1,sds2,sds3,sds4),
                   model=as.factor(rep(c(1,2,3,4), each=length(tt[-1]))),
                   weights=rep(w, each=length(tt[-1])))

# model-averaged estimates
quant <- qt(0.975, df=10)
adat <- data.frame(time=tt[-1],
                   est=apply(cbind(est1,est2,est3,est4), 1, function(x) sum(w*x)),
                   lower=apply(cbind(est1-quant*sds1,
                                     est2-quant*sds2,
                                     est3-quant*sds3,
                                     est4-quant*sds4), 1, function(x) sum(w*x)),
                   upper=apply(cbind(est1+quant*sds1,
                                     est2+quant*sds2,
                                     est3+quant*sds3,
                                     est4+quant*sds4), 1,function(x) sum(w*x)))


### graphic: time vs. delta delta c(t)
library(ggplot2)
ggplot(adat, aes(x=time, y=est)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill=grey(0.85)) +
  geom_line(size=1.3) + geom_line(data=edat, aes(group=model), linetype=2) +
  theme_bw() +
  geom_hline(yintercept=1, linetype=3) +
  xlab("Hours after resupply of phosphate to the nutrient solution") +
  ylab(expression(frac(c(t)[eEf]^(time==h) / c(t)[PT6]^(time==h), c(t)[eEf]^(time==0) / c(t)[PT6]^(time==0))))



####################################################################
####################################################################
pt6c$Time <- as.factor(pt6c$Time)
pt6result <- qpcr_nlme(response="Fluorescence", cycle="Cycle",
gene="Target", control_gene="eEf", treatment="Time",
control_treatment="0", brep="Content", well="Well", data=pt6c,
cutoff=100, conf.level=0.95, adjusted=TRUE, ratio_ddct=TRUE,
verbose=TRUE)
}

}

\keyword{datasets}
