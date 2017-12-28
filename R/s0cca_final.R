#' Final s0cca model calculation
#' @param x,y: data matrix
#' @param lambda: thresholding parameter
#' @param covI: If ture, using identity matrix for cov(x) and cov(y).
#' @return estimated coefficients and correlation
#' @export
s0cca.final=function(x,y,lambda,covI=F){

  n=dim(x)[1]
  p=dim(x)[2]
  q=dim(y)[2]
  cor_xy=cor(x,y)
  cor_xy=ifelse(abs(cor_xy)<lambda,0,cor_xy)
  cor_x=cor(x)
  cor_y=cor(y)
  if(covI){
    cov_xroot=diag(p)
    cov_yroot=diag(q)
  }else{
    if(is.positive.definite(cor_x)){
      cov_xroot=solve(sqrt.mat(cor_x))
    }else{cov_xroot=diag(p)}
    if(is.positive.definite(cor_y)){
      cov_yroot=solve(sqrt.mat(cor_y))
    }else{cov_yroot=diag(q)}
  }
  k=cov_xroot%*%cor_xy%*%cov_yroot
  k_svd=svd(k)
  u_hat=-cov_xroot%*%k_svd$u[,1]
  v_hat=-cov_yroot%*%k_svd$v[,1]
  corr=cor(x%*%u_hat,y%*%v_hat)
  return(list(u=u_hat,v=v_hat,corr=corr))
}
