## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----fig.width=6,fig.height=6-------------------------------------------------
## some unnassociated snps
set.seed(42)
nsnps=100
z1=rnorm(nsnps)
z2=rnorm(nsnps)
## spike in some strong effects
spikes=42:43
z1[spikes] = z1[spikes] + 35
z2[spikes] = z2[spikes] + 35
## add weaker effects at other snps
z1[-spikes] = z1[-spikes] + 
  rnorm(nsnps,
        mean=pmax(30-abs(1:nsnps-mean(spikes)),0),sd=1)[-spikes]
z2[-spikes] = z2[-spikes] + 
  rnorm(nsnps,
        mean=pmax(30-abs(1:nsnps-mean(spikes)),0),sd=1)[-spikes]
s1=s2=rep(sqrt(0.001),nsnps)
beta1=z1 * s1
beta2=z2 * s2
summary(exp(beta1))
summary(exp(beta2))

library(coloc)
D1=list(beta=beta1,varbeta=s1^2,type="cc",snp=paste0("v",1:100),position=1:100)
D2=list(beta=beta2,varbeta=s2^2,type="cc",snp=paste0("v",1:100),position=1:100)
oldpar=par(mfrow=c(2,1))
plot_dataset(D1); title(main="trait 1"); 
abline(v=spikes[2])
plot_dataset(D2); title(main="trait 2"); 
abline(v=spikes[2])
par(oldpar)

## -----------------------------------------------------------------------------
  coloc.abf(D1,D2)

## ----fig.width=6,fig.height=6-------------------------------------------------
variable_information=1/(sqrt(0.92))
variable_information
s1[spikes]=s1[spikes] * (c(1,variable_information))
  s2[spikes]=s2[spikes] * (c(variable_information,1))
A1=list(beta=beta1,varbeta=s1^2,type="cc",snp=paste0("v",1:100),position=1:100)
A2=list(beta=beta2,varbeta=s2^2,type="cc",snp=paste0("v",1:100),position=1:100)

oldpar <- par(mfrow = c(1,2))
plot_dataset(A1); title(main="trait 1, modified"); 
abline(v=spikes[2])
plot_dataset(A2); title(main="trait 2"); 
abline(v=spikes[2])
par(oldpar)

## -----------------------------------------------------------------------------
coloc.abf(A1,A2)

## ----fig.width=6,fig.height=6-------------------------------------------------
oldpar=par(mfrow=c(2,1))
plot(1:nsnps, finemap.abf(D1)$SNP.PP[1:nsnps], 
     xlab="Position",ylab="Posterior prob")
title(main="trait 1, original"); abline(v=spikes[2])
plot(1:nsnps, finemap.abf(D2)$SNP.PP[1:nsnps], 
     xlab="Position",ylab="Posterior prob")
title(main="trait 2, original"); abline(v=spikes[2])
par(mfrow=c(2,1))
plot(1:nsnps, finemap.abf(A1)$SNP.PP[1:nsnps], 
     xlab="Position",ylab="Posterior prob")
title(main="trait 1, modified"); abline(v=spikes[2])
plot(1:nsnps, finemap.abf(A2)$SNP.PP[1:nsnps], 
     xlab="Position",ylab="Posterior prob")
title(main="trait 2, modified"); abline(v=spikes[2])
par(oldpar)

## -----------------------------------------------------------------------------
library(colocPropTest)
LD=diag(nsnps); dimnames(LD)=list(D1$snp,D1$snp)
result=run_proptests(D1,D2,LD=LD)
min(result$fdr)

## -----------------------------------------------------------------------------
resultA=run_proptests(A1,A2,LD=LD)
min(resultA$fdr)

## ----fig.width=6,fig.height=6-------------------------------------------------
beta1=c(1,.9)
vbeta1=matrix(c(.1,.09,.09,.1),2,2)
beta2=c(1,0)
vbeta2=matrix(c(.1,0,0,.1),2,2)
library(plotrix)
plot(beta1, beta2, xlim=c(-.1,1.1), ylim=c(-.1,1.1),
     xlab=c("study 1"), ylab=c("study 2"),asp=1)
abline(0,1,lty=2)
points(0,0,pch="x"); 
text(c(beta1,0)-.05, c(beta2,0), c("snp 1","snp 2","origin"), adj=c(1,0.5))
draw.circle(beta1[1],beta2[1],.196)
draw.circle(beta1[2],beta2[2],.196); 
title(main="Effect sizes and confidence intervals")

