#' Generate Structure Matrix in gscca
#' @param edges: the edges of graph structure, see peacakge igraph
#' @param x: data matrix
#' @param plain: whether using weighted penalty. If True, the penalty is $\Sigma |\beta_i-\beta_j|$. If false, the peanlty is $\Sigma |r_ij||\beta_i-sign(r_ij)\beta_j|$.
#' @return structure matrix, D
#' @export
StrucMat=function(edges,x,plain=T){
  gr = graph(edges=edges,directed=FALSE)
  D = getDgSparse(gr)
  if(plain){return(D)}else{
    n=nrow(x);p=ncol(x);E=nrow(D)
    cor.x=cor(x)
    cor.x=cor.x+abs(min(cor.x))
    #cor.x.scaled=abs(cor.x)+1
    S=matrix(rep(0,E*p),nrow=E,ncol=p)
    for(i in 1:E){
      idx=which(D[i,]!=0)
      S[i,]=D[i,]*abs(cor.x[idx[1],idx[2]])
      #S[i,idx[1]]=sign(cor.x[idx[1],idx[2]])*S[i,idx[1]]
    }
    #S=rbind(S,diag(p)*gamma)
    S=as(S,'dgCMatrix')
    return(S)
  }
}
