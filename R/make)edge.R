#' Generate edges given tips
#' @param tip: tips, a matrix, the number of column is 2.
#' @return edges for generating structural matrix
#' @export
#'
Make_edge=function(tip){
  eg=c()
  for(i in 1:nrow(tip)){
    eg=c(eg,as.numeric(combn(tip[i,1]:tip[i,2],2)))
  }
  return(eg)
}
