#' TPR and FPR
#' @param label: true label, vector, elements are either 0 or 1
#' @param pred: predicted value, vector, elements are either 0 or 1
#' @return TPR, FPR
#' @export
TruePR=function(label,pred){
  return(sum(label*pred)/sum(label))
}
FalsePR=function(label,pred){
  return(sum((pred-label)==1)/sum(1-label))
}
