#' ROC curve, s1cca
#' @return a list containing TPR,FPR
#' @export
s1cca_roc=function(x,y,u,v,penaltyx,penaltyy){
  score=c()
  for(i in penaltyx){
    for(j in penaltyy){
      model=CCA(x,y,penaltyx = i,penaltyz = j,trace = F)
      u_hat = 1-(model$u==0);v_hat = 1-(model$v==0)
      trueu = 1-(u==0);truev=1-(v==0)
      score=rbind(score,c(TruePR(trueu,u_hat),FalsePR(trueu,u_hat),TruePR(truev,v_hat),FalsePR(truev,v_hat),
                          TruePR(c(trueu,truev),c(u_hat,v_hat)),FalsePR(c(trueu,truev),c(u_hat,v_hat))))

    }
  }
  return(list(TPRu=score[,1],FPRu=score[,2],TPRv=score[,3],FPRv=score[,4],TPR=score[,5],FPR=score[,6]))
}
