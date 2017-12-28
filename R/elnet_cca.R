#' Elastic net cca
#' @param a.u,a.v: The elasticnet mixing parameter
#' @param l.u,l.v: lambda*penalty
#' @return estimated coefficients, optimal prameters
#' @export

elnet.cca=function(x,y,a.u,a.v,l.u,l.v){
  n=nrow(x);p=ncol(x);q=ncol(y)
  cor.xy=cor(x,y)
  svdxy=svd(cor.xy)
  u0=svdxy$u[,1];v0=svdxy$v[,1]
  y.u=cor.xy%*%v0
  y.v=t(cor.xy)%*%u0
  u.hat=glmnet(diag(p),y.u,alpha = a.u,lambda = l.u)$beta
  v.hat=glmnet(diag(q),y.v,alpha = a.v,lambda = l.v)$beta
  return(list(u=u.hat,v=v.hat))
}


elnet.cca.bic=function(x,y,alpha.u,alpha.v,cv.method='ebic'){
  n=nrow(x);p=ncol(x);q=ncol(y)
  cor.xy=cor(x,y)
  svdxy=svd(cor.xy)
  u0=svdxy$u[,1];v0=svdxy$v[,1]
  y.u=cor.xy%*%v0
  y.v=t(cor.xy)%*%u0

  score=c()
  for(i in 1:length(alpha.u)){
    u.mod=glmnet(diag(p),y.u,alpha = alpha.u[i],nlambda = 20)
    lambda.u=u.mod$lambda
    for(j in 1:length(alpha.v)){
      v.mod=glmnet(diag(q),y.v,alpha=alpha.v[j],nlambda = 20)
      lambda.v=v.mod$lambda
      for(k in 1:length(lambda.u)){
        for(l in 1:length(lambda.v)){
          ri=cor(x%*%u.mod$beta[,k],y%*%v.mod$beta[,l])
          if(cv.method=='bic'){
            score.ui=bic.cca(n,u.mod$beta[,k],ri)
            score.vi=bic.cca(n,v.mod$beta[,l],ri)
          }
          if(cv.method=='ebic'){
            score.ui=ebic.cca(n,u.mod$beta[,k],ri)
            score.vi=ebic.cca(n,v.mod$beta[,l],ri)
          }
          if(cv.method=='hbic'){
            score.ui=hbic.cca(n,u.mod$beta[,k],ri)
            score.vi=hbic.cca(n,v.mod$beta[,l],ri)
          }
          if(cv.method=='gic'){
            score.ui=gic.cca(n,u.mod$beta[,k],ri)
            score.vi=gic.cca(n,v.mod$beta[,l],ri)
          }
          score=rbind(score,c(i,j,k,l,score.ui,score.vi,ri))
        }
      }
    }
  }
  i.u.opt=which.min(score[,5])
  i.v.opt=which.min(score[,6])
  a.u.opt=alpha.u[score[i.u.opt,][1]]
  a.v.opt=alpha.v[score[i.v.opt,][2]]
  l.u.opt=lambda.u[score[i.u.opt,][3]]
  l.v.opt=lambda.v[score[i.v.opt,][4]]

  u.hat=as.matrix(glmnet(diag(p),y.u,alpha = a.u.opt,lambda = l.u.opt)$beta)
  v.hat=as.matrix(glmnet(diag(q),y.v,alpha = a.v.opt,lambda = l.v.opt)$beta)
  return(list(u=u.hat,v=v.hat,corr=cor(x%*%u.hat,y%*%v.hat),
              best.param=c(a.u.opt,a.v.opt,l.u.opt,l.v.opt)))
}
