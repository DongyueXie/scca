#' Function to generate structure matrix using ape tree object
#' @description This function takes tree object from ape package as input.
#' @param x: data matrix
#' @param tree: tree object
#' @param p_c: |r_ij|^p_c
#' @param p_d: d^p_d
#' @param dis: whether distance weighted;
#' @param wcorr: whether using correlation weighted
#' @param h: controls the complexity of structural penalty. Larger h, more complex. h should between 0 and 1.
#' @export
GenStrucMat=function(x,tree,p_c=1,
                     dis=F,p_d=1,h=0.1,
                     wcorr=F,thresh=0.5){
  p=dim(x)[2]
  distMat=cophenetic(tree)
  diag(distMat)=1
  cor_x=cor(x)
  #cor_x=cor.x+abs(min(cor.x))
  cor_x=ifelse(cor_x<=thresh,0,1)
  #thresholdig the distMat to control the complexity of structural penalty
  distMat = ifelse(distMat>quantile(distMat,h),Inf,distMat)
  if(dis==F){
    distMat = ifelse(distMat!=Inf,1,Inf)
    p_d=1
  }
  #get the weight
  if(wcorr){
    weightMat=abs(cor_x)^p_c/(distMat)^p_d
  }else{
    weightMat=1/(distMat)^p_d
  }

  #generate the 'Laplacian' matrix
  comb=combn(1:p,2)
  ncomb=dim(comb)[2]
  L=matrix(rep(0,ncomb*p),nrow = ncomb,ncol = p)
  for(i in 1:ncomb){
    idx=comb[,i]
    L[i,idx[1]] = weightMat[idx[1],idx[2]]
    L[i,idx[2]] = -weightMat[idx[1],idx[2]]
  }
  L=L[which(rowSums(abs(L))!=0),]
  return(as(L,'dgCMatrix'))
}
