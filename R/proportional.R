globalVariables(c(".","snp1","snp2","p","chisquare","logp","fdr"))

##' Proportional colocalisation testing supplying only a pair of regression
##' coefficients.
##'
##' @return a list, containing
##' * result: the test statistic
##' * plot.data: dataset for plotting the input data
##' * plot.eta: dataset for plotting chisq as a function of theta or eta
##' @author Chris Wallace
##' @export
##' @param b1 regression coefficients for trait 1, expect length(b1)=2
##' @param b2 regression coefficients for trait 2, expect length(b2)=2
##' @param V1 2x2 variance-covariance matrix for trait 1
##' @param V2 2x2 variance-covariance matrix for trait 2
estprop <- function(b1,b2,V1,V2) {
    nsnps <- length(b1)
    if(nsnps!=2)
        stop("for more than 2 snps, use colocPropTest:::estprop_slow. estprop uses some results assuming exactly two snps for speed")
    stopifnot(length(b2)==nsnps,
              nrow(V1)==nsnps,
              ncol(V1)==nsnps,
              nrow(V2)==nsnps,
              ncol(V2)==nsnps)
    theta.min <- 0
    theta.max <- pi
    
    ## long winded way of writing, but allows vectorisation
    Va=matrix(c(V1[2,2],-V1[1,2],-V1[1,2],V1[1,1]),2,2)
    Wa=matrix(c(V2[2,2],-V2[1,2],-V2[1,2],V2[1,1]),2,2)
    detV=det(V1)
    detW=det(V2)
    bVb=as.vector(tcrossprod(b1,Va) %*% b1)
    bWb=as.vector(tcrossprod(b1,Wa) %*% b1)
    dVd=as.vector(tcrossprod(b2,Va) %*% b2)
    dWd=as.vector(tcrossprod(b2,Wa) %*% b2)
    bVd=as.vector(tcrossprod(b1,Va) %*% b2)
    bWd=as.vector(tcrossprod(b1,Wa) %*% b2)
    VW=sum( Wa * V1 )
    chisq_v3 <- function(theta) {
        numer=cos(theta)^4 * (dWd+bVb-bWb-dVd) +
            sin(theta)*cos(theta)^3 * 2 * (bVd-bWd) +
            cos(theta)^2 * (bWb + dVd - 2*bVb) +
            sin(theta)*cos(theta) * (-2) * bVd +
            bVb
        denom=cos(theta)^4 * ( detV + detW - VW) +
            cos(theta)^2 * ( -2*detV + VW ) +
            detV
        numer/denom
    }
    
    findmin_v3 <- function() {
        ## there are at most two minima, and never both on the same side of pi/2
        o.left <- optimize(chisq_v3,interval=c(0,pi/2))
        o.right <- optimize(chisq_v3,interval=c(pi/2,pi))
        if(o.left$objective < o.right$objective) {
            return(o.left)
        } else {
            return(o.right)
        }
    }
    
    fm <- findmin_v3()
    theta.hat <- fm$minimum; eta.hat=tan(theta.hat)
    X2 <- fm$objective
    
    x <- seq(theta.min,theta.max,length=1001)
    plot.eta <- list(theta=x,
                     eta=tan(x),
                     chisq=chisq_v3(x))
    
    list(result=c(eta.hat=eta.hat,chisquare=X2,n=nsnps,marg1=max(abs(b1/sqrt(diag(V1)))),
                  marg2=max(abs(b2/sqrt(diag(V2))))),
         plot.data=list(coef1=b1,coef2=b2,var1=diag(V1),var2=diag(V2)),
         plot.eta=plot.eta)
}

##' Proportional colocalisation testing
##'
##' This should return the same as estprop for a pair of snps, but is
##' slower. Left here for checking. Also accomodates more than two
##' snps.
##' 
##' @return a list, containing the test statistic and two datasets for
##'     plotting the input data or eta
##' @author Chris Wallace
##' @param b1 regression coefficients for trait 1
##' @param b2 regression coefficients for trait 2
##' @param V1 variance-covariance matrix for trait 1
##' @param V2 variance-covariance matrix for trait 2
estprop_slow <- function(b1,b2,V1,V2) {
    nsnps <- length(b1)
    stopifnot(length(b2)==nsnps,
              nrow(V1)==nsnps,
              ncol(V1)==nsnps,
              nrow(V2)==nsnps,
              ncol(V2)==nsnps)
    theta.min <- 0
    theta.max <- pi
    
    d <- function(theta,b1,b2) { sin(theta) * b1 - cos(theta) * b2 }
    Vstar <- function(theta) { sin(theta)^2 * V1 + cos(theta)^2 * V2 }
    chisq <- function(theta,b1,b2) { t(d(theta,b1,b2)) %*% solve(Vstar(theta)) %*% d(theta,b1,b2) }
    chisqV <- Vectorize(chisq, "theta")
    
    findmin <- function(b1,b2) {
        ## there are at most two minima, and never both on the same side of pi/2
        o.left <- optimize(chisq,interval=c(0,pi/2),b1=b1,b2=b2)
        o.right <- optimize(chisq,interval=c(pi/2,pi),b1=b1,b2=b2)
        if(o.left$objective < o.right$objective) {
            return(o.left)
        } else {
            return(o.right)
        }
    }
    
    fm <- findmin(b1,b2)
    theta.hat <- fm$minimum; eta.hat=tan(theta.hat)
    X2 <- fm$objective
    
    x <- seq(theta.min,theta.max,length=1001)
    plot.eta <- list(theta=x,
                     eta=tan(x),
                     chisq=chisqV(x,b1,b2))
    
    list(result=c(eta.hat=eta.hat,chisquare=X2,n=nsnps,marg1=max(abs(b1/sqrt(diag(V1)))),
                  marg2=max(abs(b2/sqrt(diag(V2))))),
         plot.data=list(coef1=b1,coef2=b2,var1=diag(V1),var2=diag(V2)),
         plot.eta=plot.eta)
}

#' create variance-covariance matrix for pair of marginal beta + vbeta, given
#' estimate of r between snps
#'
#' @param beta vector of two coefficients at two snps
#' @param vbeta vector of two variance of coefficients at the same two snps
#' @param rho LD (r) between the two snps
#'
#' @return list of coefficient & variance-covariance matrix
marg_with_V=function(beta,vbeta,rho) {
    res=list(beta=beta)
    res$V=diag(vbeta)
    res$V[1,2]=res$V[2,1]=prod(sqrt(diag(res$V))) * rho
    res
}

#' run proportional test directly on marginal test stats from coloc datasets
#'
#' @param j indices of thw two snps
#' @param S1 coloc dataset 1
#' @param S2 coloc dataset 2
#' @param LD1 LD matrix for dataset 1 - rownames, colnames capture the snps and S1$snp[j] must
#'   be represented
#' @param LD2 LD matrix for dataset 2 - rownames, colnames capture the snps and S2$snp[j] must
#'   be represented. if not supplied, defaults to LD1
#'
#' @return result from estprop
tester_marg=function(j,S1,S2,LD1,LD2=LD1) {
    k=j[2]
    j=j[1]
    snpj=S1$snp[ j ]
    snpk=S1$snp[ k ]
    j1=marg_with_V(beta=S1$beta[c(j,k)],
                   vbeta=S1$varbeta[c(j,k)],
                   rho=LD1[snpj,snpk])
    j2=marg_with_V(beta=S2$beta[c(j,k)],
                   vbeta=S2$varbeta[c(j,k)],
                   rho=LD2[snpj,snpk])
    res=estprop(as.vector(j1$beta), as.vector(j2$beta), j1$V, j2$V)
    res$result
}

#' Helper function to adjust LD parameter r for differential sample size between snps
#' 
#' Estimate the r between effect estimates at snps which were
#' genotyped on different sets of cases and controls.  The adjusted r
#' will be nform(...) * r (where r is the population correlation
#' between snps).
#'
#' @param n0a number of controls with data at snp a
#' @param n1a number of cases with data at snp a
#' @param n0b number of controls with data at snp b
#' @param n1b number of cases with data at snp b
#' @param n0ab number of controls with data at both snps a and b
#' @param n1ab number of cases with data at both snps a and b
#'
#' @return proportionality constant that depends on sample size.
#'
nform=function(n0a,n1a,n0b,n1b,n0ab=pmin(n0a,n0b),n1ab=pmin(n1a,n1b)) {
    f=sqrt(n1a * n1b / n0a / n0b)
    (n1ab / f + n0ab * f) / sqrt(n1a+n0a) / sqrt(n1b+n0b)
}


#' adjust LD for variable sample size
#'
#' @param S coloc style dataset, with additional entries n0  and n1  which are *vectors* giving the number of cases and controls genotyped at each SNP
#' @param LD matrix of LD, with dimnames given by snps in S$snp
#' @return adjusted LD matrix
#' @export
#'
#' @examples
#' library(coloc)
#' 
#' data(coloc_test_data)
#' attach(coloc_test_data)
#' LD=D1$LD
#' dimnames(LD)=list(D1$snp,D1$snp)
#' D1$type="cc"
#' D1$s=.5
#' D1$n1=D1$N * sample(c(0.25,.5),length(D1$snp), replace=TRUE)
#' D1$n0=rep(0.5*D1$N,length(D1$snp))
#' aLD=colocPropTest::adjust_LD(D1,LD)
#' LD[1:6,1:6]
#' aLD[1:6,1:6]
#' detach(coloc_test_data)
adjust_LD=function(S,LD) {
    ## check
    check_dataset(S)
    stopifnot("n0, n1 in S"="n1" %in% names(S) && "n0" %in% names(S),
              "n0, n1 vectors with one entry per snp"=length(S$n1)==length(S$snp) && length(S$n0)==length(S$snp))
    snp1=rownames(LD)[ row(LD) ]
    snp2=colnames(LD)[ col(LD) ]
    n1a=S$n1[ match(snp1, S$snp) ]
    n0a=S$n0[ match(snp1, S$snp) ]
    n1b=S$n1[ match(snp2, S$snp) ]
    n0b=S$n0[ match(snp2, S$snp) ]
    adj=nform(n0a=n0a,n0b=n0b,n1a=n1a,n1b=n1b)
    LD * matrix(adj, nrow(LD), ncol(LD))
}

##' run proportional tests on extreme subset of snp pairs from two
##' coloc style datasets. Of all functions in this package, this is
##' the main one that should be used.
##'
##' @param S1 coloc dataset 1
##' @param S2 coloc dataset 2
##' @param LD LD matrix - rownames, colnames capture the snps and S1$snp[j] must
##'   be represented
##' @param topsnps list of topsnps to be considered for testing or, if "auto",
##'   will be automatically selected
##' @param r2.thr r2 threshold for initial tagging step - includes only one of
##'   any set of snps in mutually high LD with r2 > r2.thr
##' @param maxtests maximum number of test pairs to consider. if more than
##'   maxtests pairs available, will select a random sample.
##' @param nauto number of snps to use when automatically defining topsnps. only
##'   has an effect if topsnps=="auto"
##' @param adjust_n TRUE if you want to adjust for variable sample size between
##'   snps. This is only set up for case control data at the moment (ask if you
##'   need quantitative) and requires that you supply separately the number of
##'   cases and controls at each snp in each dataset, as vector elements of the lists
##'   called n1 (cases) and n0 (controls)
##' @return data.table containing the tests run
##' @export
##' @author Chris Wallace
##'
##' @examples
##' library(colocPropTest)
##' library(coloc)
##' data(coloc_test_data)
##' attach(coloc_test_data)
##' LD=D1$LD
##' dimnames(LD)=list(D1$snp,D1$snp)
##' results=run_proptests(D1,D2,LD=LD,topsnps=D1$snp,maxtests=100)
##' min(results$fdr)
run_proptests=function(S1,S2,LD,topsnps="auto",r2.thr=0.95,maxtests=10000, nauto=200, adjust_n=FALSE) {
    if(length(topsnps)==1 && topsnps=="auto") {
        z1=abs(S1$beta/sqrt(S1$varbeta))
        z2=abs(S2$beta/sqrt(S2$varbeta))
        snps_both=intersect(S1$snp, S2$snp)
        z_both=log(z1[ match(snps_both, S1$snp) ]) + log(z2[ match(snps_both, S2$snp) ])
        topsnps=head(snps_both[ order(z_both, decreasing=TRUE) ], nauto)
    }
    ## LD pruning, using modified tag() from chr1swallace/GUESSFM
    tags=tag(LD[topsnps,topsnps]^2,r2.thr)
    message("after tagging, unique topsnps: ",length(tags))
                                        # pairs is index of snp in topsnps
    pairs=expand.grid(which(topsnps %in% tags), which(topsnps %in% tags)) %>% t()
    pairs=pairs[, pairs[1,]<pairs[2,] ]
                                        # pairs is now character using snp ids
                                        #pairs=rbind( topsnps[ pairs[1,] ],
                                        #             topsnps[ pairs[2,] ])
    dim(pairs) # 3828
    
    n=length(intersect(S1$snp,S2$snp)) # this is used to define how many tests we could have done
    LD=LD[ topsnps, topsnps ] 
    ok=LD[t(pairs)] < sqrt(0.99) # don't do tests where LD unknown or snps are in rsq > r2.thr
    pairs=pairs[,ok]
    if(ncol(pairs) > maxtests) {
        use=sample(1:ncol(pairs), maxtests)
        pairs=pairs[,use]
    }
    if(adjust_n) {
        LD1=adjust_LD(LD, S1)
        LD2=adjust_LD(LD, S2)
    } else {
        LD1=LD2=LD
    }
    message("possible tests: ",n*(n-1)/2)
    message("tests to run: ",ncol(pairs)) 
    if(ncol(pairs)<2) {
        warning("no strongly significant tests at pairs of snps in incomplete LD")
        return(NULL)
    }
    
    S1=keep_from_S(S1,topsnps)
    S2=keep_from_S(S2,topsnps)
    
    stats=lapply(1:ncol(pairs), function(i) {
                                        #message(i)
        tester_marg(pairs[,i],S1=S1,S2=S2,LD1=LD1,LD2=LD2)
    })
    
    result = stats %>% do.call("rbind",.) %>% as.data.frame() %>% as.data.table()
                                        # this isn't needed now we are directly using the marginal stats - it caught issues with failure to estimate joint models. but it doesn't do harm, so leave it just in case
    nulls=sapply(stats,is.null)
    if(any(nulls))
        warning("some null results in proportional testing: ",sum(nulls), " / ", length(nulls))
    result[,snp1:=topsnps[ pairs[1,] ][ !nulls ]]
    result[,snp2:=topsnps[ pairs[2,] ][ !nulls ]]
    result[,p:=pchisq(chisquare,1,lower.tail=FALSE)]
    result[,logp:=-pchisq(chisquare,1,lower.tail=FALSE,log.p=TRUE)]
    ntests=n*(n-1)/2
    result[,fdr:=p.adjust(p, n=ntests)]
    attr(result,"ntests")=ntests
    attr(result,"nrun")=length(stats)
    attr(result,"nnulls")=sum(nulls)
    return(result)
}
