#' Calculate Matrix^(-1/2)
#' @param m: input matrix
#' @return m^(-1/2)
#' @export
sqrt.mat=function(m){
  m=nearPD(m,corr=T)$mat
  eigen.m=eigen(m)
  ev=(eigen.m$values)^(-1/2)
  ev[is.na(ev)]=0
  return(eigen.m$vectors%*%diag(ev)%*%t(eigen.m$vectors))
}
