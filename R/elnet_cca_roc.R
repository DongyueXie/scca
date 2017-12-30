#' ROC curve, elnet
#' @return a list containing TPR,FPR
#' @export
elnet.cca.roc=function(x,y,u,v,alpha.u,alpha.v,nlambda=20){
  n=nrow(x);p=ncol(x);q=ncol(y)
  cor.xy=cor(x,y)
  svdxy=svd(cor.xy)
  u0=svdxy$u[,1];v0=svdxy$v[,1]
  y.u=cor.xy%*%v0
  y.v=t(cor.xy)%*%u0

  score=c()
  for(i in 1:length(alpha.u)){
    u.mod=glmnet(diag(p),y.u,alpha = alpha.u[i],nlambda = nlambda)
    lambda.u=u.mod$lambda
    for(j in 1:length(alpha.v)){
      v.mod=glmnet(diag(q),y.v,alpha=alpha.v[j],nlambda = nlambda)
      lambda.v=v.mod$lambda
      for(k in 1:length(lambda.u)){
        for(l in 1:length(lambda.v)){
          u_hat = 1-(u.mod$beta[,k]==0);v_hat = 1-(v.mod$beta[,l]==0)
          trueu = 1-(u==0);truev=1-(v==0)
          score=rbind(score,c(TruePR(trueu,u_hat),FalsePR(trueu,u_hat),TruePR(truev,v_hat),FalsePR(truev,v_hat),
                              TruePR(c(trueu,truev),c(u_hat,v_hat)),FalsePR(c(trueu,truev),c(u_hat,v_hat))))

        }
      }
    }
  }
  return(list(TPRu=score[,1],FPRu=score[,2],TPRv=score[,3],FPRv=score[,4],TPR=score[,5],FPR=score[,6]))
}
