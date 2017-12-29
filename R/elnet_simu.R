#' Simulation for elastic net cca
#' @param alpha.u,alpha.v: The elasticnet mixing parameter, between 0 and 1.
#' @export
elnet.cca.simu=function(niter,u,v,n=50,data.type='normal',
                        sigma_z=1,sigma_x=NULL,sigma_y=NULL,
                        theta=0.01,n_x=1000,n_y=1000,lambda_z=10,
                        cv.method='bic',alpha.u=c(0.6,0.8,1),alpha.v=c(0.6,0.8,1),
                        K=5,verbose=F){

  p=length(u);q=length(v)
  corrs=c()
  us=c()
  vs=c()
  para_opt=c()
  TPR_u=c()
  TPR_v=c()
  FPR_v=c()
  FPR_u=c()
  TPR=c()
  FPR=c()
  F1_u=c()
  F1_v=c()
  F1=c()
  for(iter in 1:niter){
    #print(iter)
    set.seed(iter+1)
    if(data.type=='normal'){
      datax = GenStrucData_Normal(u=u,v=v,n=n,sigma_z = sigma_z,sigma_x = sigma_x,sigma_y = sigma_y)
    }
    if(data.type=='count_normal'){
      datax = GenStrucData_Count_normal(u=u,v=v,n=n,sigma_z = sigma_z,
                                        sigma_x = sigma_x,sigma_y = sigma_y,theta=theta,n_x=n_x,n_y=n_y)
    }
    if(data.type=='count_latent'){
      datax = GenStrucData_Count_latent(u=u,v=v,n=n,theta=theta,lambda_z = lambda_z,
                                        n_x=n_x,n_y=n_y)
    }
    x=datax$x
    y=datax$y
    if(cv.method=='bic'){
      elnet.model=elnet.cca.bic(x,y,alpha.u = alpha.u,alpha.v = alpha.v,cv.method = cv.method)
    }
    if(cv.method=='ebic'){
      elnet.model=elnet.cca.bic(x,y,alpha.u = alpha.u,alpha.v = alpha.v,cv.method = cv.method)
    }
    if(cv.method=='hbic'){
      elnet.model=elnet.cca.bic(x,y,alpha.u = alpha.u,alpha.v = alpha.v,cv.method = cv.method)
    }
    if(cv.method=='gic'){
      elnet.model=elnet.cca.bic(x,y,alpha.u = alpha.u,alpha.v = alpha.v,cv.method = cv.method)
    }
    if(cv.method=='cv'){
      elnet.model=elnet.cca.cv(x,y,alpha.u = alpha.u,alpha.v = alpha.v,K=K)
    }


    para_opt=cbind(para_opt,elnet.model$best.param)
    corrs[iter]=elnet.model$corr
    u_hat=elnet.model$u;v_hat=elnet.model$v
    #record u and v
    us=cbind(us,u_hat)
    vs=cbind(vs,v_hat)
    u_hat = 1-(u_hat==0);v_hat = 1-(v_hat==0)
    trueu = 1-(u==0);truev=1-(v==0)
    TPR_u[iter]=sum(u_hat*trueu)/sum(trueu)
    TPR_v[iter]=sum(v_hat*truev)/sum(truev)
    FPR_u[iter]=sum((u_hat-trueu)==1)/sum(1-trueu)
    FPR_v[iter]=sum((v_hat-truev)==1)/sum(1-truev)
    TPR[iter]=(sum(u_hat*trueu)+sum(v_hat*truev))/(sum(trueu)+sum(truev))
    FPR[iter]=(sum((u_hat-trueu)==1)+sum((v_hat-truev)==1))/(sum(1-trueu)+sum(1-truev))
    F1_u[iter]=2*sum(u_hat*trueu)/(2*sum(u_hat*trueu)+sum((u_hat-trueu)==1)+sum((u_hat-trueu)==-1))
    F1_v[iter]=2*sum(v_hat*truev)/(2*sum(v_hat*truev)+sum((v_hat-truev)==1)+sum((v_hat-truev)==-1))
    F1[iter]=(2*sum(u_hat*trueu)+2*sum(v_hat*truev))/((2*sum(u_hat*trueu)+sum((u_hat-trueu)==1)+sum((u_hat-trueu)==-1))+(2*sum(v_hat*truev)+sum((v_hat-truev)==1)+sum((v_hat-truev)==-1)))
    if(verbose&iter%%10==0){print(sprintf('Iteration %s, TPR is %.2f, FPR is %.2f.',iter,TPR[iter],FPR[iter]))}
  }
  return(list(n=n,p=p,q=q,corr=corrs,us=us,vs=vs,para_opt=para_opt,
              TPR=TPR,TPR_u=TPR_u,TPR_v=TPR_v,
              FPR=FPR,FPR_u=FPR_u,FPR_v=FPR_v,F1_u=F1_u,F1_v=F1_v,F1=F1))

}
