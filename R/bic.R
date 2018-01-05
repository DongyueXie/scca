#' BIC functions for CCA
#'
#' @param n: sample size
#' @param b: estimated coefficients, vector
#' @param r: canonical correlation coefficient
#' @param g: a parameter in gbic and hbic
#' @return bic, a scalar
#' @export
bic.cca=function(n,b,r){
  return(n*log(1-r^2)+sum(b!=0)*log(n))
}

ebic.cca=function(n,b,r,g=0.5){
  return(n*log(1-r^2)+sum(b!=0)*log(n)+2*sum(b!=0)*g*log(length(b)))
}

hbic.cca=function(n,b,r,g=1.5){
  return(n*log(1-r^2)+2*sum(b!=0)*g*log(length(b)))
}

gic.cca=function(n,b,r){
  return(n*log(1-r^2)+sum(b!=0)*length(b)^(1/3))
}
