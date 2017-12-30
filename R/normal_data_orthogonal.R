#' Generate normal data, orthogonal direction
#'
#' @return data x and y
#' @export
#'
GenStrucData_Normal_oth=function(u,v,n=50,sigma_z=1,sigma_x=NULL,sigma_y=NULL){
  p=length(u)
  q=length(v)
  if(is.null(sigma_x)){sigma_x=diag(1,p,p)}
  if(is.null(sigma_y)){sigma_y=diag(1,q,q)}
  x=c()
  y=c()
  u1=u;v1=v
  u2=0.5*(-u1);v2=0.5*(-v1)
  for(i in 1:n){
    zi1=rnorm(1,0,sigma_z)
    zi2=rnorm(1,0,sigma_z)
    xi=zi1*u1+zi2*u2+mvrnorm(1,rep(0,p),sigma_x)
    x=rbind(x,xi)
    yi=zi1*v1+zi2*v2+mvrnorm(1,rep(0,q),sigma_y)
    y=rbind(y,yi)
  }
  #print (sigma_x, sigma_y, sigma_z)
  return(list(x=x,y=y))
}
