#' ROC curve, s0cca
#' @return a list containing TPR,FPR
#' @export

s0cca.roc = function(x,y,u,v,ncluster=NULL,cluster.type='no',maxsteps=20){
  #library(matrixcalc)
  #library(expm)
  n=dim(x)[1]
  p=dim(x)[2]
  q=dim(y)[2]
  covxy=cor(x,y)
  covxy_abs=abs(covxy)
  cc=apply(abs(covxy),2,max)
  rr=apply(abs(covxy),1,max)
  lambda=sort(unique(c(cc,rr)),decreasing = T)
  #
  if(is.null(ncluster)){ncluster=length(lambda)}
  #how to choose the lambdas:seq,kmeans,quantile
  if(cluster.type=='seq'){
    lambda=seq(max(lambda),min(lambda),length.out = ncluster)
  }
  if(cluster.type=='quantile'){
    choose.lambda=c()
    for(i in 1:ncluster){
      choose.lambda[i]=quantile(lambda,(ncluster+1-i)/ncluster)
    }
    lambda=choose.lambda
  }
  if(cluster.type=='kmeans'){
    cc=kmeans(lambda,ncluster-1)
    lambda=sort(cc$centers,decreasing = T)
  }
  lambda=lambda[-1]

  #maxstrp controls how many lambdas to search from.
  if(is.null(maxsteps)){maxsteps=length(lambda)}

  score=c()

    #the other cv methods
  for(i in 1:maxsteps){
      #print(i)
      covxy_hat=ifelse(covxy_abs<=lambda[i],0,covxy)
      u.proxy=rowSums(abs(covxy_hat))
      v.proxy=colSums(abs(covxy_hat))
      u_hat = 1-(u.proxy==0);v_hat = 1-(v.proxy==0)
      trueu = 1-(u==0);truev=1-(v==0)
      score=rbind(score,c(TruePR(trueu,u_hat),FalsePR(trueu,u_hat),TruePR(truev,v_hat),FalsePR(truev,v_hat),
                          TruePR(c(trueu,truev),c(u_hat,v_hat)),FalsePR(c(trueu,truev),c(u_hat,v_hat))))
  }
  return(list(TPRu=score[,1],FPRu=score[,2],TPRv=score[,3],FPRv=score[,4],TPR=score[,5],FPR=score[,6]))
}

