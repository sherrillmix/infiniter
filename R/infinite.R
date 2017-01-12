#' Generate fake cluster data
#' @param nCluster number of clusters
#' @param nPoints number of points within each cluster
#' @param nDim number of dimensions
#' @param betweenSd cluster centers are normally distributed around 0 with this standard deviation
#' @param withinSd points are normally distributed around cluster centers with this standard deviation
#' @return list with matrix giving point locations (rows are points, columns are dimensions), a matrix of cluster centers and vector of cluster IDs
#' @author Scott Sherrill-Mix \email{R@@sherrillmix.com}
#' @export
#' @examples
#' clusts<-generateFakeClusters()
#' plot(clusts$points,col=clusts$ids)
#' points(clusts$centers,pch='x',col=1:nrow(clusts$centers),cex=2)
generateFakeClusters<-function(nCluster=2,nPoints=100,nDim=2,betweenSd=3,withinSd=1){
  centers<-lapply(1:nCluster,function(xx)stats::rnorm(nDim,0,betweenSd))
  points<-mapply(function(center,n)do.call(rbind,replicate(n,stats::rnorm(length(center),center,withinSd),simplify=FALSE)),centers,nPoints,SIMPLIFY=FALSE)
  return(list('points'=do.call(rbind,points),'centers'=do.call(rbind,centers),'ids'=rep(1:nCluster,sapply(points,nrow))))
}


infiniteGaussian<-function(x,alpha=1,initK=3,nIter=10){
  n<-nrow(x)
  z<-sample(initK,n)
  K<-initK
  for(iter in 1:nIter){
    zTab<-table(z)
    for(ii in 1:n){
      thisTab<-zTab
      thisTab[as.character(z[ii])]<-thisTab[as.character(z[ii])]-1
      pZi<-thisTab/(n+alpha-1)
      pXi<-
      for(kk in 1:K){
      } 
    }
  }
}
