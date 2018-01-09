#' ROC curve, gscca
#' @return a list containing TPR,FPR
#' @export

gscca.roc=function(x,y,u,v,edgex,edgey,maxsteps=20,plain=T,
                   gamma.u,gamma.v,lambda.u=NULL,lambda.v=NULL,Sx=NULL,Sy=NULL){
  n=dim(x)[1]
  init=gscca(x=x,y=y,edgex=edgex,edgey=edgey,maxsteps=maxsteps,plain=plain,Sx=Sx,Sy=Sy)
  out.u=init$out.u
  out.v=init$out.v
  if(is.null(lambda.u)){lambda.u=out.u$lambda}
  if(is.null(lambda.v)){lambda.v=out.v$lambda}
  score=c()
  for(i in 1:length(gamma.u)){
    beta.ui=softthresh(out.u,lambda.u,gamma.u[i])

    for(j in 1:length(gamma.v)){
      beta.vi=softthresh(out.v,lambda.v,gamma.v[j])

      for(k in 1:dim(beta.ui)[2]){
        for(l in 1:dim(beta.vi)[2]){
          u_hat = 1-(beta.ui[,k]==0);v_hat = 1-(beta.vi[,l]==0)
          trueu = 1-(u==0);truev=1-(v==0)
          score=rbind(score,c(TruePR(trueu,u_hat),FalsePR(trueu,u_hat),TruePR(truev,v_hat),FalsePR(truev,v_hat),
                              TruePR(c(trueu,truev),c(u_hat,v_hat)),FalsePR(c(trueu,truev),c(u_hat,v_hat))))
        }
      }
    }
  }
  return(list(TPRu=score[,1],FPRu=score[,2],TPRv=score[,3],FPRv=score[,4],TPR=score[,5],FPR=score[,6]))
}
