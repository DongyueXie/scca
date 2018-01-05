#' Return cca correlation
#' @param x,y: data
#' @param ratio: n/p, when to use identity matrix as cor
#' @export

cancor2=function(x,y,ratio=2){
  n=nrow(x)
  p=ncol(x)
  q=ncol(y)
  if(n/p>=ratio){
    Mat1=sqrt.mat(cor(x))
  }else{Mat1=diag(p)}
  if(n/q>=ratio){
    Mat2=sqrt.mat(cor(y))
  }else{Mat2=diag(q)}
  K=svd(Mat1%*%cor(x,y)%*%Mat2)
  return(cor(x%*%K$u[,1],y%*%K$v[,1]))
}
