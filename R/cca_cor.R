#' Re-estimate the canonical correlaiton
#' @param x,y: data
#' @param u,v: estiamted coefficients
#' @export

cca.cor=function(x,y,u,v){
  m1=try(cancor2(x[,u!=0],y[,v!=0]),T)
  if(class(m1)=='try-error'){
    return(cor(x%*%u,y%*%v))
  }else{
    return(m1)
  }
}
