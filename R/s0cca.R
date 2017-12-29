#' Conduct s0cca
#' @param x,y: data matrix
#' @param ncluster: the number of cluster when clustering  the lambda. If NULL, using all lambdas.
#' @param cluster.type: how to cluster lambdas, choose from 'no'(do not do the clustering), 'seq','quantile','kmeans'.
#' @param maxsteps: how many lambdas to try, default to be 20.
#' @param cv.method: choose from 'cv','bic','ebic','hbic','gic'.
#' @return optimal canonical coeffcients, canonical correlation, optimal lambda.
#' @export


s0cca = function(x,y,ncluster=NULL,cluster.type='no',maxsteps=20,cv.method='cv'){
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

  cv.u=c()
  cv.v=c()
  us=c()
  vs=c()
  rs=c()

  idx=KFold(1:n,5)

  if(cv.method=='cv'){
    corr_vall=c()
    for(i in 1:maxsteps){
      #print(i)
      covxy_hat=ifelse(covxy_abs<=lambda[i],0,covxy)
      u.proxy=rowSums(abs(covxy_hat))
      v.proxy=colSums(abs(covxy_hat))
      corr_val=c()
      for(j in 1:length(idx)){
        #print(j)
        x_train=x[-idx[[j]],]
        x_val = x[idx[[j]],]
        y_train=y[-idx[[j]],]
        y_val = y[idx[[j]],]
        train=s0cca.thresh(x_train,y_train,u.proxy,v.proxy)
        x_val=x_val[,u.proxy!=0]
        y_val=y_val[,v.proxy!=0]
        x_val=cbind(x_val,rep(0,length(idx[[j]])))
        y_val=cbind(y_val,rep(0,length(idx[[j]])))
        corr_val[j]=cor(x_val%*%train$u,y_val%*%train$v)

      }
      corr_vall[i]=mean(corr_val)
    }
    ll=lambda[which.max(corr_vall)]
    finalmodel=s0cca.final(x,y,ll)
    return(list(u=finalmodel$u,v=finalmodel$v,corr=finalmodel$corr,lambda=ll))
  }else{
    #the other cv methods
    for(i in 1:maxsteps){
      #print(i)
      covxy_hat=ifelse(covxy_abs<=lambda[i],0,covxy)
      u.proxy=rowSums(abs(covxy_hat))
      v.proxy=colSums(abs(covxy_hat))
      us=cbind(us,u.proxy)
      vs=cbind(vs,v.proxy)
      covxy_hat=covxy_hat[u.proxy!=0,v.proxy!=0]

      svdi=svd(covxy_hat)
      ui=svdi$u[,1]
      vi=svdi$v[,1]
      if(sum(u.proxy!=0)==1&sum(v.proxy!=0)==1){
        r=cor(x[,u.proxy!=0]*ui,y[,v.proxy!=0]*vi)
      }else if(sum(u.proxy!=0)==1){
        svdi=svd(t(covxy_hat))
        ui=svdi$u[,1]
        vi=svdi$v[,1]
        cor.y=cor(y[,v.proxy!=0])
        is.pos.y=is.positive.definite(cor.y)
        if(is.pos.y){
          r=cor(x[,u.proxy!=0]*ui,y[,v.proxy!=0]%*%sqrt.mat(cor.y)%*%vi)
        }else{
          r=cor(x[,u.proxy!=0]*ui,y[,v.proxy!=0]%*%vi)
        }

      }else if(sum(v.proxy!=0)==1){

        cor.x=cor(x[,u.proxy!=0])
        is.pos.x=is.positive.definite(cor.x)
        if(is.pos.x){
          r=cor(x[,u.proxy!=0]%*%sqrt.mat(cor.x)%*%ui,y[,v.proxy!=0]*vi)
        }else{
          r=cor(x[,u.proxy!=0]%*%ui,y[,v.proxy!=0]*vi)
        }

      }else{
        cor.x=cor(x[,u.proxy!=0])
        cor.y=cor(y[,v.proxy!=0])
        is.pos.x=is.positive.definite(cor.x)
        is.pos.y=is.positive.definite(cor.y)
        if(is.pos.x&is.pos.y){
          r=cor(x[,u.proxy!=0]%*%sqrt.mat(cor.x)%*%ui,y[,v.proxy!=0]%*%sqrt.mat(cor.y)%*%vi)
        }else if(is.pos.x){
          r=cor(x[,u.proxy!=0]%*%sqrt.mat(cor.x)%*%ui,y[,v.proxy!=0]%*%vi)
        }else if(is.pos.y){
          r=cor(x[,u.proxy!=0]%*%ui,y[,v.proxy!=0]%*%sqrt.mat(cor.y)%*%vi)
        }else{
          r=cor(x[,u.proxy!=0]%*%ui,y[,v.proxy!=0]%*%vi)
        }
      }

      rs[i]=r
      if(cv.method=='bic'){
        cv.u[i]=bic.cca(n,ui,r)
        cv.v[i]=bic.cca(n,vi,r)
      }
      if(cv.method=='ebic'){
        cv.u[i]=ebic.cca(n,ui,r)
        cv.v[i]=ebic.cca(n,vi,r)
      }
      if(cv.method=='hbic'){
        cv.u[i]=hbic.cca(n,ui,r)
        cv.v[i]=hbic.cca(n,vi,r)
      }
      if(cv.method=='gic'){
        cv.u[i]=gic.cca(n,ui,r)
        cv.v[i]=gic.cca(n,vi,r)
      }

    }
    i.u=which.min(cv.u)
    i.v=which.min(cv.v)
    i.r=which.max(rs)
    if(cv.method=='bic'|cv.method=='hbic'|cv.method=='ebic'|cv.method=='gic'){
      return(list(u=us[,i.u],v=vs[,i.v],corr=rs[i.r]))
    }else{
      return(list(u=us[,i.r],v=vs[,i.r],corr=rs[i.r]))
    }
  }

}
