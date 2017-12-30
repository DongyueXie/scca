#' Moving average
#'
#' @param x: data
#' @param n: average interval
#' @return averaged data
#' @export
#'
moav=function(x,n){return(filter(x,rep(1/n,n),sides=2))}
