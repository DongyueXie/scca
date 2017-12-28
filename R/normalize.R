#' Normalize a random variable
#'
#' @param x: input data matrix
#' @return normalized column, matrix
#' @export
normalize = function(x){
  return(apply(x,2,function(y){(y-mean(y))/sd(y)}))
}
