#' Generate count data
#' @description Generate count data, using normal data from latent variable model
#' @param u: caonical coefficients for x
#' @param v: canonical coefficients for y
#' @param n: sample size
#' @param sigma_z: sd of the latent variable z, scalar
#' @param sigma_x: covariance matrix of x, default to be identity
#' @param sigma_y: covariance matrix of y, default to be identity
#' @param n_y: parameter n in DM model
#' @param n_x: parameter n in DM model
#' @param theta: oversidpersion parameter in DM model
#' @return a list containing data matrix X and Y
#' @export

GenStrucData_Count_normal = function(u,v,n=50,n_y=1000,n_x=1000,theta=0.01,
                                     sigma_z=1,sigma_x=NULL,sigma_y=NULL){
  datax=GenStrucData_Normal(u,v,n,sigma_z,sigma_x,sigma_y)
  x=datax$x;y=datax$y
  x=apply(x,1,function(x){x1=(x-min(x))/(max(x)-min(x));return(x1/sum(x1))})
  y=apply(y,1,function(x){x1=(x-min(x))/(max(x)-min(x));return(x1/sum(x1))})
  x1=apply(x,1,function(x){return(simPop(1,n=n_x,pi=x,theta=theta)$data)})
  y1=apply(y,1,function(x){return(simPop(1,n=n_y,pi=x,theta=theta)$data)})
  return(list(x=x1,y=y1))
}
