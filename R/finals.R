#' Summary of simulation results
#' @param object: an object from simulation function
#' @return TPR,FPR,F1,Cor
#' @export
#'
finals=function(object){
  results=round(cbind(
    mean(object$TPR_u),
    sd(object$TPR_u),
    mean(object$FPR_u),
    sd(object$FPR_u),
    mean(object$TPR_v),
    sd(object$TPR_v),
    mean(object$FPR_v),
    sd(object$FPR_v),
    mean(object$TPR),
    sd(object$TPR),
    mean(object$FPR),
    sd(object$FPR),
    mean(object$F1_u),
    sd(object$F1_u),
    mean(object$F1_v),
    sd(object$F1_v),
    mean(object$F1),
    sd(object$F1),
    mean(object$corr)),4)
  colnames(results)=c('TPR_u','sd_TPR_u','FPR_u','sd_FPR_u','TPR_v','sd_TPR_v','FPR_v','sd_FPR_v',
                      'TPR','sd_TPR','FPR','sd_FPR','F1_u','sd_F1_u','F1_v','sd_F1_v','F1','sd_F1','Cor')
  return(results)
}
