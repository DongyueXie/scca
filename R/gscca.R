#' Implement gscca model
#' @param minlam: the minimum lambda at which path algo should stop.
#' @return canonical correlation coefficients
#' @export
gscca=function(x,y,edgex,edgey,maxsteps=2000,minlam.u=NULL,minlam.v=NULL,plain=T,Sx=NULL,Sy=NULL,thresh=0.5){
  #initialize
  n=nrow(x);p=ncol(x);q=ncol(y)
  cor.xy=cor(x,y)
  svdxy=svd(cor.xy)
  u0=svdxy$u[,1];v0=svdxy$v[,1]
  y.u=cor.xy%*%v0
  y.v=t(cor.xy)%*%u0
  if(is.null(Sx)){
    Sx=StrucMat(edgex,x,plain,thresh)
    if(identical(fusedlasso(y.u,D=Sx,maxsteps=2)$lambda,0)){
      Sx=StrucMat(edgex,x,T,thresh)
    }
  }else{
    if(identical(fusedlasso(y.u,D=Sx,maxsteps=2)$lambda,0)){
      message('Sx entries are too small')
    }
  }
  if(is.null(Sy)){
    Sy=StrucMat(edgey,y,plain,thresh)
    if(identical(fusedlasso(y.v,D=Sy,maxsteps=2)$lambda,0)){
      Sy=StrucMat(edgey,y,T,thresh)
    }
  }else{
    if(identical(fusedlasso(y.v,D=Sy,maxsteps=2)$lambda,0)){
      message('Sy entries are too small')
    }
  }
  #
  if(is.null(minlam.u)){
    out.u=fusedlasso(y.u,D=Sx,maxsteps=maxsteps)
  }else{
    out.u=fusedlasso(y.u,D=Sx,minlam = minlam.u)
  }
  if(is.null(minlam.v)){
    out.v=fusedlasso(y.v,D=Sy,maxsteps=maxsteps)
  }else{
    out.v=fusedlasso(y.v,D=Sy,minlam = minlam.v)
  }

  #
  return(list(out.u=out.u,out.v=out.v))

}
