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

#' Fit an infinite Bayesian mixture model to IPD and sequence data
#' @param ipds a matrix of IPD ratios with a row for each position of interest and a column for each surrounding base
#' @param seqs a matrix of single characters or vector of strings giving sequences of interest
#' @param nIter number of JAGS iterations
#' @param maxK the maximum number of clusters expected (smaller is less computationally expensive)
#' @param mc.cores number of cores (and chains) to use
#' @return list with matrix giving simulated data and stored JAGS model (can be used for further iterations)
#' @author Scott Sherrill-Mix \email{R@@sherrillmix.com}
#' @export
#' @examples
#' 1
infiniteGaussian<-function(ipds,seqs,nIter=100,maxK=100,mc.cores=4){
  n<-nrow(ipds)
  if(nrow(seqs)!=n)stop('Numbers of seqs does not match ipds')
  if(!is.matrix(seqs)){
    if(any(nchar(seqs)!=nchar(seqs[1])))stop('seqs not all same length')
    seqMat<-do.call(rbind,strsplit(seqs,''))
  }else{
    seqMat<-seqs
  }
  jagsModel <- '
    var pwm[nSeq,4,k];
    model {
      a ~ dunif(0.3, 100)
      pwmAlpha ~ dunif(.3,100)
      for (ii in 1:k){
        alpha[ii] <- a/k
        for(jj in 1:nIpd){
          mu[ii,jj] ~ dnorm(0.0, 1)
          sigma[ii,jj] ~ dunif(0, 10)
        }
        for(jj in 1:nSeq){
          pwm[jj,,ii] ~ ddirch(c(pwmAlpha,pwmAlpha,pwmAlpha,pwmAlpha))
        }
      }
      p[1:k] ~ ddirich(alpha[1:k])
      for (ii in 1:n){
        z[ii] ~ dcat(p)
        for(jj in 1:nIpd){
          ipd[ii,jj] ~ dnorm(mu[z[ii],jj], sigma[z[ii],jj])
        }
        for(jj in 1:nSeq){
          seq[jj,,ii] ~ dmulti(pwm[jj,,z[ii]],1)
        }
      }
    }
  '
  chains<-parallel::mclapply(sample(1:1e6,mc.cores),function(xx){
    model <- rjags::jags.model(
    textConnection(jagsModel),
    data=list(ipd=ipds, n=n, k=maxK,nIpd=ncol(ipds),nSeq=ncol(seqMat),seq=seqMat),
    inits=list(.RNG.seed=as.integer(xx),'.RNG.name'='base::Wichmann-Hill'),
    n.adapt=nIter
    )
    chain <-rjags::coda.samples(model = jagsModel, n.iter = nIter, variable.names = c('p', 'mu', 'sigma','z','pwm','pwmAlpha','a'))
    rchain <- as.matrix(chain)
    return(list('sims'=rchain,'model'=model))
  },mc.cores=mc.cores)
}
