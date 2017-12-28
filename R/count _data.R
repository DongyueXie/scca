#' Generate count data
#' @description Generate count data using latent variable model
#' @param lambda_z: poisson parameter of the latent variable z
#' @param theta: oversidpersion parameter in DM model
#' @param n_y: parameter n in DM model
#' @param n_x: parameter n in DM model
#' @param u: caonical coefficients for x
#' @param v: canonical coefficients for y
#' @param n: sample size
#' @export
GenStrucData_Count_latent = function(u,v,n=50,lambda_z=10,theta=0.01,
                                     n_x=1000,n_y=1000){

  p=length(u)
  q=length(v)

  x=c()
  y=c()
  for(i in 1:n){
    zi=rpois(1,lambda_z)
    pxi = rep(1/p,p)
    xi=zi*u+simPop(1,pi=pxi,n=n_x,theta=theta)$data
    x=rbind(x,xi)
    pyi = rep(1/q,q)
    yi=zi*v+simPop(1,pi=pyi,n=n_y,theta=theta)$data
    y=rbind(y,yi)
  }
  return(list(x=x,y=y))
}
