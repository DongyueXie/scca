#' Permutation Test
#' @description Perform permutaiton test of gscca
#' @param nperm: the number of permutation test
#' @return p.corr,p.corr.re: p values
#' @return corrs,corrs.re: all correlations from permutation
#' @return corr0,corr.re0: optimal correlation
#' @export
#'
gscca.permute=function(nperm=100,x,y,edgex,edgey,maxsteps=20,plain=T,cv.method='bic',
                       gamma.u,gamma.v,lambda.u,lambda.v,
                       Sx=NULL,Sy=NULL,ccor=F,thresh=0.5,seed=NULL){
  set.seed(seed)
  n=dim(x)[1]
  p=dim(x)[2]
  q=dim(y)[2]
  corrs = c()
  corrs.re=c()
  for(pm in 1:nperm){
    permute = sample(1:n,n,replace = F)
    x_p = x[permute,]
    mod.fit=gscca.bic(x_p,y,edgex,edgey,maxsteps,plain,cv.method,
                      gamma.u,gamma.v,lambda.u,lambda.v,
                      Sx,Sy,ccor,thresh)
    corrs[pm]=ifelse(is.na(mod.fit$corr),0,mod.fit$corr)
    corrs.re[pm]=ifelse(is.na(mod.fit$corr.re),0,mod.fit$corr.re)
  }
  mod0=gscca.bic(x,y,edgex,edgey,maxsteps,plain,cv.method,
                 gamma.u,gamma.v,lambda.u,lambda.v,
                 Sx,Sy,ccor,thresh)
  corr0=ifelse(is.na(mod0$corr),0,mod0$corr)
  corr.re0=ifelse(is.na(mod0$corr.re),0,mod0$corr.re)
  p.corr=mean(c(corrs)>=c(corr0))
  p.corr.re=sum(c(corrs.re)>=c(corr.re0))/nperm
  return(list(p.corr=p.corr,p.corr.re=p.corr.re,corrs=corrs,corrs.re=corrs.re,
              corr0=corr0,corr.re0=corr.re0))
}
