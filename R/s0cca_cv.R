#' CV function for s0cca
#' @description For each fold of CV, return the estimated coefficients and correlation
#' @param x,y: data
#' @param u,v: vector indcating the sparseness of coefficients
#' @return estimated coefficients and correlation
#' @export
s0cca.thresh=function(x,y,u,v){
  n=dim(x)[1]
  p=dim(x)[2]
  q=dim(y)[2]
  x.spar=x[,u!=0]
  y.spar=y[,v!=0]
  x.spar=cbind(x.spar,rep(0,n))
  y.spar=cbind(y.spar,rep(0,n))
  cor.x=cor(x.spar)
  cor.x[is.na(cor.x)]=0
  cor.y=cor(y.spar)
  cor.y[is.na(cor.y)]=0
  cor.xy=cor(x.spar,y.spar)
  cor.xy[is.na(cor.xy)]=0
  cor.x.r=sqrt.mat(cor.x)
  cor.y.r=sqrt.mat(cor.y)
  #print(cor.x.r);print(cor.y.r)
  k_svd=try(svd(cor.x.r%*%cor.xy%*%cor.y.r),TRUE)
  if(class(k_svd)=='try-error'){k_svd=svd(cor.xy)}
  u_hat=cor.x.r%*%k_svd$u[,1]
  v_hat=cor.y.r%*%k_svd$v[,1]
  corr=k_svd$d[1]
  return(list(u=u_hat,v=v_hat,corr=corr))
}
