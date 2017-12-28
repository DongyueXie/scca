#' gscca model selection using bic
#' @param gamma,lambda: lambda*fused penalty+gamma*lambda*l1 penalty.
#' @return estimated coefficients, correlation and optimal parameters
#' @export
gscca.bic=function(x,y,edges,maxsteps=2000,plain=T,cv.method='bic',
                   gamma.u,gamma.v,lambda.u=NULL,lambda.v=NULL,Sx=NULL,Sy=NULL){
  n=dim(x)[1]
  init=gscca(x=x,y=y,edges=edges,maxsteps=maxsteps,plain=plain,Sx=Sx,Sy=Sy)
  out.u=init$out.u
  out.v=init$out.v
  if(is.null(lambda.u)){lambda.u=out.u$lambda}
  if(is.null(lambda.v)){lambda.v=out.v$lambda}
  #all the u and v estimated without 1 norm penalty.
  #beta.u=coef(out.u,lambda.u)$beta
  #beta.v=coef(out.v,lambda.v)$beta
  #cross validation
  score=c()
  for(i in 1:length(gamma.u)){
    beta.ui=softthresh(out.u,lambda.u,gamma.u[i])
    #the sparsity penalty meybe too large
    if(sum(beta.ui)==0){beta.ui=beta.ui+0.001}
    for(j in 1:length(gamma.v)){
      beta.vi=softthresh(out.v,lambda.v,gamma.v[j])
      if(sum(beta.vi)==0){beta.vi=beta.vi+0.001}
      for(k in 1:dim(beta.ui)[2]){
        for(l in 1:dim(beta.vi)[2]){
          ri=cor(x%*%beta.ui[,k],y%*%beta.vi[,l])
          if(cv.method=='bic'){
            score.ui=bic.cca(n,beta.ui[,k],ri)
            score.vi=bic.cca(n,beta.vi[,l],ri)
          }
          if(cv.method=='ebic'){
            score.ui=ebic.cca(n,beta.ui[,k],ri)
            score.vi=ebic.cca(n,beta.vi[,l],ri)
          }
          if(cv.method=='hbic'){
            score.ui=hbic.cca(n,beta.ui[,k],ri)
            score.vi=hbic.cca(n,beta.vi[,l],ri)
          }
          if(cv.method=='gic'){
            score.ui=gic.cca(n,beta.ui[,k],ri)
            score.vi=gic.cca(n,beta.vi[,l],ri)
          }
          score=rbind(score,c(i,j,k,l,score.ui,score.vi,ri))
        }
      }
    }
  }
  lambda.u.opt=lambda.u[score[which.min(score[,5]),][3]]
  lambda.v.opt=lambda.v[score[which.min(score[,6]),][4]]
  gamma.u.opt=gamma.u[score[which.min(score[,5]),][1]]
  gamma.v.opt=gamma.v[score[which.min(score[,6]),][2]]
  #print(score)
  #print(score[which.min(score[,5]),][1])
  #print(score[which.min(score[,6]),][2])
  u.hat=softthresh(out.u,lambda.u.opt,gamma.u.opt)
  v.hat=softthresh(out.v,lambda.v.opt,gamma.v.opt)
  return(list(u=u.hat,v=v.hat,corr=cor(x%*%u.hat,y%*%v.hat),
              best.param=c(lambda.u.opt,lambda.v.opt,gamma.u.opt,gamma.v.opt)))
}

