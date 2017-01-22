#' Generate fake cluster data
#' @param nCluster number of clusters
#' @param nPoints number of points within each cluster
#' @param nIpd number of IPDs to generate
#' @param nBase number of bases to generate
#' @param betweenSd cluster centers are normally distributed around 0 with this standard deviation
#' @param withinSd points are normally distributed around cluster centers with this standard deviation
#' @param nVariantBases number of bases in each PWM varying from equal proportions
#' @return list with matrix giving point locations (rows are points, columns are dimensions), a matrix of cluster centers, a vector of sequences, a list of PWMs (rows are positions, columns are bases) and vector of cluster IDs
#' @author Scott Sherrill-Mix \email{R@@sherrillmix.com}
#' @export
#' @examples
#' clusts<-generateFakeClusters()
#' plot(clusts$ipds,col=clusts$ids)
#' points(clusts$centers,pch='x',col=1:nrow(clusts$centers),cex=2)
generateFakeClusters<-function(nCluster=2,nPoints=100,nIpd=6,nBase=3,nVariantBases=min(3,nBase),betweenSd=3,withinSd=1){
  centers<-lapply(1:nCluster,function(xx)stats::rnorm(nIpd,0,betweenSd))
  ipds<-mapply(function(center,n)do.call(rbind,replicate(n,stats::rnorm(length(center),center,withinSd),simplify=FALSE)),centers,nPoints,SIMPLIFY=FALSE)
  pwms<-mapply(function(x,nVar){
    pwm<-matrix(.25,ncol=4,nrow=nBase)
    varBases<-sample(nBase,nVar)
    for(ii in varBases){
      pwm[ii,]<-sample(rStickBreak(4))
    }
    return(pwm)
  },1:nCluster,nVariantBases,SIMPLIFY=FALSE)
  bases<-c('A','C','G','T')
  seqs<-apply(do.call(rbind,mapply(function(pwm,n){
      do.call(cbind,lapply(1:nBase,function(ii)sample(bases,n,TRUE,prob=pwm[ii,])))
  },pwms,nPoints,SIMPLIFY=FALSE)),1,paste,collapse='')
  return(list('ipds'=do.call(rbind,ipds),'centers'=do.call(rbind,centers),'seqs'=seqs,'pwms'=pwms,'ids'=rep(1:nCluster,sapply(ipds,nrow))))
}

#' Convert sequences to a numbered matrix
#' @param seqs a character vector of sequences
#' @param bases the potential characters in the sequences in the order to be numbered
#' @return 3-d logical array with the first dimension being sequence, second position and third base
#' @author Scott Sherrill-Mix \email{R@@sherrillmix.com}
#' @export
#' @examples
#' seqToMat(c('AAAT','ATCG'))
seqToMat<-function(seqs,bases=c('A','C','G','T')){
  if(any(nchar(seqs)!=nchar(seqs[1])))stop('All sequences not same length')
  seqMat<-do.call(rbind,strsplit(seqs,''))
  return(do.call(abind::abind,c(lapply(bases,function(x)seqMat==x),list(along=3))))
}


#' Calculate random proportions from a dirchlet process based on stick breaking procedure
#'
#' @param n number of numbers to generate
#' @param alpha smaller alpha means on average more weight concentrated earlier
#' @return a vector of n numbers adding to 1 
#' @export
#' @examples
#' rStickBreak(10)
#' rStickBreak(10,3)
rStickBreak<-function(n,alpha=1){
  sticks<-c()
  for(ii in 1:(n-1))sticks<-c(sticks,stats::rbeta(1,1,alpha)*(1-sum(sticks)))
  sticks<-c(sticks,1-sum(sticks))
  return(sticks)
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
#' fakeSeqs<-generateFakeClusters()
#' infiniteGaussian(fakeSeqs$ipds,fakeSeqs$seqs,maxK=10)
infiniteGaussian<-function(ipds,seqs,nIter=100,maxK=100,mc.cores=4){
  n<-nrow(ipds)
  if(length(seqs)!=n)stop('Numbers of seqs does not match ipds')
  seqMat<-seqToMat(seqs)
  jagsModel <- '
    var pwm[k,nBase,4];
    model {
      a ~ dunif(0.3, 100)
      pwmAlpha ~ dunif(.3,100)
      for (ii in 1:k){
        alpha[ii] <- a/k
        for(jj in 1:nIpd){
          mu[ii,jj] ~ dnorm(0.0, 1)
          sigma[ii,jj] ~ dunif(0, 10)
        }
        for(jj in 1:nBase){
          pwm[ii,jj,] ~ ddirch(c(pwmAlpha,pwmAlpha,pwmAlpha,pwmAlpha))
        }
      }
      p[1:k] ~ ddirich(alpha[1:k])
      for (ii in 1:n){
        z[ii] ~ dcat(p)
        for(jj in 1:nIpd){
          ipd[ii,jj] ~ dnorm(mu[z[ii],jj], sigma[z[ii],jj])
        }
        for(jj in 1:nBase){
          seq[ii,jj,] ~ dmulti(pwm[z[ii],jj,],1)
        }
      }
    }
  '
  chains<-parallel::mclapply(sample(1:1e6,mc.cores),function(xx){
    model <- rjags::jags.model(
    textConnection(jagsModel),
    data=list(ipd=ipds, n=n, k=maxK,nIpd=ncol(ipds),nBase=ncol(seqMat),seq=seqMat),
    inits=list(.RNG.seed=as.integer(xx),'.RNG.name'='base::Wichmann-Hill'),
    n.adapt=nIter
    )
    chain <-rjags::coda.samples(model = model, n.iter = nIter, variable.names = c('p', 'mu', 'sigma','z','pwm','pwmAlpha','a'))
    rchain <- as.matrix(chain)
    return(list('sims'=rchain,'model'=model))
  },mc.cores=mc.cores)
  return(chains)
}

#' Run additional sampling on infininte Gaussian model
#'
#' @param chains the output from infiniteGaussian
#' @param nIter number of additional iterations
#' @export
#' @return list with matrix giving simulated data and stored JAGS model (can be used for further iterations)
iterateIG<-function(chains,nIter=500){
  out<-parallel::mclapply(chains,function(xx){
    sims<-rjags::coda.samples(model = xx, n.iter = nIter, variable.names = c('p', 'mu', 'sigma','z','pwm','pwmAlpha','a'))
    return(list('sims'=sims,'model'=chains))
  })
  return(out)
}
