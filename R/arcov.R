#' Generate autoregressive cov matrix
#' @param r: correlation for neighboring variable
#' @param p: dimension
#' @return a cov matrix
#' @export
arcov=function(r,p){
  sigma_x = matrix(rep(0,p^2),p,p)
  for(i in 1:p){
    for(j in 1:p){
      sigma_x[i,j]=(r)^abs(i-j)
    }
  }
  return(sigma_x)
}
