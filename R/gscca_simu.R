#' Simulation for gscca
#' @param gamma.u,gamma.v: ratio of l1 peanlty and fused penalty
#' @param lambda.u,lambda.v: fused penalty. lambda*fused penalty+gamma*lambda*l1 penalty.
#' @param Sx,Sy: structure matrix for Sx,Sy.
#' @param maxsteps: the number of lambdas to try, default to be 20. Largest 20 values from solution path
#' @param ccor: whether re-estiamte the coeeficient when doing BIC
#' @export
gscca.simu.normal=function(niter,edgex,edgey,u,v,n=50,data.type='normal',
                           sigma_z=1,sigma_x=NULL,sigma_y=NULL,
                           theta=0.01,lambda_z=10,n_x=1000,n_y=1000,
                           cv.method='bic',gamma.u=NULL,gamma.v=NULL,
                           lambda.u=NULL,lambda.v=NULL,maxsteps=20,plain=T,
                           Sx=NULL,Sy=NULL,verbose=F,K=5,ccor=F){
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
  for(j in 1:niter){
    #print(j)
    set.seed(j+1)
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
    x=datax$x
    y=datax$y
    #use the l1 penalty parameter from s1cca
    if(is.null(gamma.u)){
      cc=CCA.permute(x,y,penaltyxs = seq(0.01,0.7,length.out = 10),
                     penaltyzs = seq(0.01,0.7,length.out = 10),trace = F)
      gamma.u=unique(c(cc$bestpenaltyx*sqrt(n),0.3*sqrt(n),cc$bestpenaltyx))
      gamma.v=unique(c(cc$bestpenaltyz*sqrt(n),0.3*sqrt(n),cc$bestpenaltyz))

      #print(gamma.u);print(gamma.v)
    }
    gscca.model=gscca.bic(x,y,edgex=edgex,edgey=edgey,gamma.u=gamma.u,gamma.v=gamma.v,
                            lambda.u = lambda.u,lambda.v = lambda.v,
                            maxsteps = maxsteps,plain = plain,Sx=Sx,Sy=Sy,cv.method=cv.method,ccor=ccor)

    # if(cv.method=='ebic'){
    #   gscca.model=gscca.bic(x,y,edges=edges,gamma.u=gamma.u,gamma.v=gamma.v,
    #                         lambda.u = lambda.u,lambda.v = lambda.v,
    #                         maxsteps = maxsteps,plain = plain,Sx=Sx,Sy=Sy,cv.method=cv.method)
    # }
    # if(cv.method=='hbic'){
    #   gscca.model=gscca.bic(x,y,edges=edges,gamma.u=gamma.u,gamma.v=gamma.v,
    #                         lambda.u = lambda.u,lambda.v = lambda.v,
    #                         maxsteps = maxsteps,plain = plain,Sx=Sx,Sy=Sy,cv.method=cv.method)
    # }
    # if(cv.method=='gic'){
    #   gscca.model=gscca.bic(x,y,edges=edges,gamma.u=gamma.u,gamma.v=gamma.v,
    #                         lambda.u = lambda.u,lambda.v = lambda.v,
    #                         maxsteps = maxsteps,plain = plain,Sx=Sx,Sy=Sy,cv.method=cv.method)
    # }

    para_opt=cbind(para_opt,gscca.model$best.param)
    corrs[j]=gscca.model$corr
    u_hat=gscca.model$u;v_hat=gscca.model$v
    #record u and v
    us=cbind(us,u_hat)
    vs=cbind(vs,v_hat)
    u_hat = 1-(u_hat==0);v_hat = 1-(v_hat==0)
    trueu = 1-(u==0);truev=1-(v==0)
    TPR_u[j]=sum(u_hat*trueu)/sum(trueu)
    TPR_v[j]=sum(v_hat*truev)/sum(truev)
    FPR_u[j]=sum((u_hat-trueu)==1)/sum(1-trueu)
    FPR_v[j]=sum((v_hat-truev)==1)/sum(1-truev)
    TPR[j]=(sum(u_hat*trueu)+sum(v_hat*truev))/(sum(trueu)+sum(truev))
    FPR[j]=(sum((u_hat-trueu)==1)+sum((v_hat-truev)==1))/(sum(1-trueu)+sum(1-truev))
    F1_u[j]=2*sum(u_hat*trueu)/(2*sum(u_hat*trueu)+sum((u_hat-trueu)==1)+sum((u_hat-trueu)==-1))
    F1_v[j]=2*sum(v_hat*truev)/(2*sum(v_hat*truev)+sum((v_hat-truev)==1)+sum((v_hat-truev)==-1))
    F1[j]=(2*sum(u_hat*trueu)+2*sum(v_hat*truev))/((2*sum(u_hat*trueu)+sum((u_hat-trueu)==1)+sum((u_hat-trueu)==-1))+(2*sum(v_hat*truev)+sum((v_hat-truev)==1)+sum((v_hat-truev)==-1)))
    if(verbose&j%%10==0){print(sprintf('Iteration %s, TPR is %.2f, FPR is %.2f.',j,TPR[j],FPR[j]))}
  }
  return(list(n=n,p=p,q=q,corr=corrs,us=us,vs=vs,para_opt=para_opt,
              TPR=TPR,TPR_u=TPR_u,TPR_v=TPR_v,
              FPR=FPR,FPR_u=FPR_u,FPR_v=FPR_v,F1_u=F1_u,F1_v=F1_v,F1=F1))
}
