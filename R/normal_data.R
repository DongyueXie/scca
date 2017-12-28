#' Generate normal data using latent variable model
#'
#' @param u: caonical coefficients for x
#' @param v: canonical coefficients for y
#' @param n: sample size
#' @param sigma_z: sd of the latent variable z, scalar
#' @param sigma_x: covariance matrix of x, default to be identity
#' @param sigma_y: covariance matrix of y, default to be identity
#' @return a list containing data mattrix X and Y
#' @export
GenStrucData_Normal=function(u,v,n=50,sigma_z=1,sigma_x=NULL,sigma_y=NULL){
  p=length(u)
  q=length(v)
  if(is.null(sigma_x)){sigma_x=diag(1,p,p)}
  if(is.null(sigma_y)){sigma_y=diag(1,q,q)}
  x=c()
  y=c()
  for(i in 1:n){
    zi=rnorm(1,0,sigma_z)
    xi=zi*u+mvrnorm(1,rep(0,p),sigma_x)
    x=rbind(x,xi)
    yi=zi*v+mvrnorm(1,rep(0,q),sigma_y)
    y=rbind(y,yi)
  }
  #print (sigma_x, sigma_y, sigma_z)
  return(list(x=x,y=y))
}
