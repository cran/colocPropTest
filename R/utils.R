##' Derive tag SNPs using heirarchical clustering
##'
##' Uses complete linkage and the \code{\link{hclust}} function to define clusters,
##' then cuts the tree at 1-tag.threshold
##'
##' @param r2 matrix of rsquared values
##' @param quiet if FALSE (default), show progress messages
##' @param r2_threshold only 1 of a set of snps with r2 > r2_threshold will be kept
##' @param method method used for heirarchical clustering.  See hclust for options.
##'
##' @return character vector, names are \code{snps}, values are the tag for each SNP
##' @author Chris Wallace
##' @export
tag <- function(r2,r2_threshold=0.95, quiet=FALSE,method="complete") {
    stopifnot("r2 should be a square matrix" = is.matrix(r2) && nrow(r2)==ncol(r2),
              "0 < r2_threshold < 1"=r2_threshold > 0 && r2_threshold<1)
    D <- as.dist(1-r2)
    hc <- hclust(D, method=method)
    clusters <- cutree(hc, h=1-r2_threshold)
    names(clusters)[!duplicated(clusters)]
}

##' draw two ellipses
##'
##' @param b1 ellipse 1 centre (2d)
##' @param vb1 ellipse 1 vcov matrix
##' @param b2 ellipse 2 centre (2d)
##' @param vb2 ellipse 2 vcov matrix
##' @param legend character vector length 2 naming ellipse 1 and 2
##' @param include_origin if TRUE, ensure plot includes (0,0)
##' @param ... arguments passed to plot()
##' @return draw ellipses on current graphics device
##' @examples
##' plot_ellipses(b1=c(5,5), vb1=diag(2),
##'               b2=c(2,2), vb2=matrix( c(1,0.5,0.5,1), 2, 2 ),
##'               legend=c("circle", "ellipse"),
##'               include.origin=TRUE) 
##' @export
##' @author Chris Wallace
plot_ellipses=function(b1,vb1,b2,vb2,legend=c("inferred","observed"),
                       include_origin=FALSE,...) {
    stopifnot("b1, b2 are length 2 vectors"=length(b1)==2 && length(b2)==2,
              "vb1, vb2 are 2x2 matrices"=is.matrix(vb1) && is.matrix(vb2) &&
                  ncol(vb1)==2 && ncol(vb2)==2 & nrow(vb1)==2 && nrow(vb2)==2)
    sets=rbind(ellipse(as.vector(b1), shape=vb1, radius=1.96,draw=FALSE),
               ellipse(as.vector(b2), shape=vb2, radius=1.96,col=1,draw=FALSE))
    mins=apply(sets,2,min)
    maxs=apply(sets,2,max)
    if(include_origin) {
        mins[ mins > 0 ] = -0.01
        maxs[ maxs < 0 ] = 0.01
    }
    plot(x=c(b1[1],b2[1]),
         y=c(b1[2],b2[2]),
         col=c("blue","black"),
         xlim=c(mins[1],maxs[1]),
         ylim=c(mins[2],maxs[2]),
         ...)
    ellipse(as.vector(b1), shape=vb1, radius=1.96)
    ellipse(as.vector(b2), shape=vb2, radius=1.96,col=1)
    legend("topleft",col=c("blue","black"),pch=c(1,1),legend=legend)
}

#' logp
#' 
#' uses logs in calculation to avoid numerical issues with very small std errors / p values
#'
#' @param beta coefficient
#' @param se std error of coefficient
#'
#' @return -log10 p
#' @export
lp=function(beta, se) {
    -(pnorm(-abs(beta/se),log.p=TRUE) + log(2))/log(10)
}

#' keep snp subset of coloc dataset
#'
#' @param S coloc dataset
#' @param keep snps to keep
#'
#' @return subset of coloc dataset
keep_from_S=function(S, keep) {
    w=match(keep,S$snp)
    if(all(is.na(w))) {
        warning("no snps matched, not subsetting")
        return(S)
    }
    if(any(is.na(w)))
        w=setdiff(w,NA)
    for(nm in which(sapply(S,is.vector) & sapply(S,length)>1))
        S[[ nm ]]  = S[[ nm ]][w]
    if("LD" %in% names(S))
        S$LD=S$LD[ w, w, drop=FALSE]
    S
}

