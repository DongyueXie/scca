#' Simulation function for s0cca
#' @param niter: the number of simulation iterations
#' @param data.type: 'normal', 'count_normal','count_latent'
#' @return TPR,FPR,F1 score
#' @export
s0cca_simu = function(niter,u,v,n=50,data.type='normal',
                      sigma_z=1,sigma_x=NULL,sigma_y=NULL,
                      theta=0.01,lambda_z=10,n_x=1000,n_y=1000,
                      ncluster=NULL,cluster.type='no',maxsteps=20,cv.method='cv',K=5,verbose=F){
  p=length(u);q=length(v)
  corrs=c()
  us=c()
  vs=c()
  lambda=c()
  TPR_u=c()
  TPR_v=c()
  FPR_v=c()
  FPR_u=c()
  F1_u=c()
  F1_v=c()
  F1=c()
  TPR=c()
  FPR=c()
  for(iter in 1:niter){
    set.seed(iter+1)
    if(data.type=='normal'){
      datax = GenStrucData_Normal(u=u,v=v,n=n,sigma_z = sigma_z,sigma_x = sigma_x,sigma_y = sigma_y)
    }
    if(data.type=='normal_oth'){
      datax = GenStrucData_Normal_oth(u=u,v=v,n=n,sigma_z = sigma_z,sigma_x = sigma_x,sigma_y = sigma_y)
    }
    if(data.type=='normal_mis'){
      datax = GenStrucData_Normal_mis(u=u,v=v,n=n,sigma_z = sigma_z,sigma_x = sigma_x,sigma_y = sigma_y)
    }
    if(data.type=='count_normal'){
      datax = GenStrucData_Count_normal(u=u,v=v,n=n,sigma_z = sigma_z,
                                        sigma_x = sigma_x,sigma_y = sigma_y,theta=theta,n_x=n_x,n_y=n_y)
    }
    if(data.type=='count_latent'){
      datax = GenStrucData_Count_latent(u=u,v=v,n=n,theta=theta,lambda_z = lambda_z,
                                        n_x=n_x,n_y=n_y)
    }

    #print(iter)
    result=s0cca(datax$x,datax$y,ncluster=ncluster,cluster.type=cluster.type,
                 maxsteps=maxsteps,cv.method=cv.method,K=K)
    corrs[iter]=result$corr
    us=cbind(result$u,us)
    vs=cbind(result$v,vs)
    lambda[iter]=result$lambda
    u_hat = 1-(result$u==0);v_hat = 1-(result$v==0)
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
  return(list(n=n,p=p,q=q,corr=corrs,us=us,vs=vs,lambda=lambda,
              TPR=TPR,TPR_u=TPR_u,TPR_v=TPR_v,
              FPR=FPR,FPR_u=FPR_u,FPR_v=FPR_v,F1_u=F1_u,F1_v=F1_v,F1=F1))
}
