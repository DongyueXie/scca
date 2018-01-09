#' Simulate the ROC
#' @description A function to run the simulation for plotting ROC
#' @param nroc: the number of iterations
#' @param data.type: type of data to generate, choose from 'normal','normal_oth','normal_oth'...
#' @param  maxstep_s0,maxstep_gs: paramters for s0cca and gscca
#' @param penaltyx,penaltyy: for s1cca, vectors
#' @param alpha.u,alpha.v,nlambda: for enscca, vector, vector, scalar.
#' @param gamma.u,gamma.v: paramters for gscca, vectors.
#' @export
roc.simu=function(nroc=1,u,v,n=50,data.type='normal',
                  sigma_z=1,sigma_x=NULL,sigma_y=NULL,
                  theta=0.01,n_x=1000,n_y=1000,lambda_z=10,
                  maxstep_s0=NULL,
                  penaltyx=seq(0,1,length.out = 50),penaltyy=seq(0,1,length.out = 50),
                  alpha.u=seq(0,1,length.out = 10),alpha.v=seq(0,1,length.out = 10),nlambda=20,
                  edgex,edgey,plain=T,maxstep_gs=20,
                  gamma.u=seq(0,5,length.out = 10),gamma.v=seq(0,5,length.out = 10),
                  seed=12345){
  mats0=c()
  mats1=c()
  matelnet=c()
  matgscca=c()
  nroc=1
  set.seed(seed)
  for(roc in 1:nroc){

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
    s0=matrix(unlist(s0cca.roc(x,y,u,v,maxsteps=maxstep_s0)),ncol=6,byrow=F)
    s1=matrix(unlist(s1cca_roc(x,y,u,v,penaltyx,penaltyy)),ncol=6,byrow=F)
    en=matrix(unlist(elnet.cca.roc(x,y,u,v,alpha.u,alpha.v,nlambda)),ncol=6,byrow=F)
    gs=matrix(unlist(gscca.roc(x,y,u,v,edgex,edgey,maxstep_gs,gamma.u = gamma.u,gamma.v = gamma.v,plain=plain)),ncol=6,byrow=F)
    mats0=rbind(mats0,s0)
    mats1=rbind(mats1,s1)
    matelnet=rbind(matelnet,en)
    matgscca=rbind(matgscca,gs)
  }
  if(nroc>1){
    mats0=apply(mats0,2,mean)
    mats1=apply(mats1,2,mean)
    matelnet=apply(matelnet,2,mean)
    matgscca=apply(matgscca,2,mean)
  }
  return(list(mats0=mats0,mats1=mats1,matelnet=matelnet,matgscca=matgscca))
}
