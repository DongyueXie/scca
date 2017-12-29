#'KFold CV Index
#' @param x: data matrix
#' @param k: number of folds
#' @return A list of cv index
#' @export
KFold=function(idxs,k){
  n=length(idxs)
  idx=list()
  #idxs=1:n
  s=floor(n/k)
  rm=n%%k
  for(i in 1:k){
    if(i<=k-rm){
      f=sample(idxs,s,replace = F)
      idx[[i]]=f
    }else{
      f=sample(idxs,s+1,replace = F)
      idx[[i]]=f
    }
    idxs=setdiff(idxs,f)
  }
  return(idx)
}
