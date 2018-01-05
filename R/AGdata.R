#' Association study for American Gut project
#' @param nOTU: the number of OTUs to use
#' @param S1,S2: choose from 'FECAL','ORAL','SKIN'
#' @param pat:path to the otu_table.biom and 97_otus.tree
#' @export
#'
AGdata=function(nOTU=100,S1,S2,pat){
  AG = read_hdf5_biom(paste(pat,'/otu_table.biom',sep = ''))
  AG_tree = read_tree_greengenes(paste(pat,'/97_otus.tree',sep = ''))
  AG_OTU = AG$data
  AG_Taxa = AG$rows
  AGmeta=readLines(paste(pat,'/ag-cleaned.txt',sep = ''))
  ag=lapply(AGmeta,function(x){unlist(strsplit(x, split = "\t"))})
  ag=data.frame(matrix(unlist(ag), nrow=length(ag), byrow=T),stringsAsFactors=FALSE)
  colnames(ag)=ag[1,]
  ag=ag[-1,]
  AGmeta=ag
  #find the paricipants that provided more than 1 samples.
  IDunique=unique(AGmeta$HOST_SUBJECT_ID)
  ID2=unique(AGmeta$HOST_SUBJECT_ID[-match(IDunique,AGmeta$HOST_SUBJECT_ID)])
  AGmeta2=AGmeta[AGmeta$HOST_SUBJECT_ID%in%ID2,]
  AGmeta2_fecals=AGmeta2[AGmeta2$SIMPLE_BODY_SITE==S1,]
  AGmeta2_orals=AGmeta2[AGmeta2$SIMPLE_BODY_SITE==S2,]
  Sub_fecal_oral=intersect(AGmeta2_fecals$HOST_SUBJECT_ID,AGmeta2_orals$HOST_SUBJECT_ID)
  AGmeta2_fecal=c()
  AGmeta2_oral=c()
  for(i in Sub_fecal_oral){
    n_fecal=sum(AGmeta2_fecals$HOST_SUBJECT_ID==i)
    n_oral=sum(AGmeta2_orals$HOST_SUBJECT_ID==i)
    n=min(c(n_fecal,n_oral))
    if(n_fecal>1|n_oral>1){
      proxy1=rbind(AGmeta2_fecals[AGmeta2_fecals$HOST_SUBJECT_ID==i,],rep(0,ncol(AGmeta2_fecals)))
      proxy2=rbind(AGmeta2_orals[AGmeta2_orals$HOST_SUBJECT_ID==i,],rep(0,ncol(AGmeta2_orals)))
      AGmeta2_fecal=rbind(AGmeta2_fecal,proxy1[1:n,])
      AGmeta2_oral=rbind(AGmeta2_oral,proxy2[1:n,])
    }else{
      AGmeta2_fecal=rbind(AGmeta2_fecal,AGmeta2_fecals[AGmeta2_fecals$HOST_SUBJECT_ID==i,])
      AGmeta2_oral=rbind(AGmeta2_oral,AGmeta2_orals[AGmeta2_orals$HOST_SUBJECT_ID==i,])
    }
  }
  ID_fecal=AGmeta2_fecal$`#SampleID`
  ID_oral=AGmeta2_oral$`#SampleID`


  ID_fecal_otu=unlist(lapply(names(AG_OTU[[1]]), function(x) x %in% ID_fecal))
  ID_oral_otu=unlist(lapply(names(AG_OTU[[1]]), function(x) x %in% ID_oral))

  #total counts for each otu
  OTUSum_fecal = sapply(AG_OTU, function(x)sum(x[ID_fecal_otu]))
  OTUSum_oral = sapply(AG_OTU, function(x)sum(x[ID_oral_otu] ))

  #top 100 counts otu's index
  OTU_fecal=OTUSum_fecal > sort(OTUSum_fecal,decreasing=TRUE)[nOTU+1]
  OTU_oral=OTUSum_oral > sort(OTUSum_oral,decreasing=TRUE)[nOTU+1]

  X=c()
  Y=c()
  for(i in 1:nrow(AGmeta2_oral)){
    fecaldata=AG_OTU[OTU_fecal]
    oraldata=AG_OTU[OTU_oral]
    sampleid_fecal=AGmeta2_fecal$`#SampleID`[i]
    sampleid_oral=AGmeta2_oral$`#SampleID`[i]
    X=rbind(X,as.numeric(sapply(fecaldata,function(x){idx=(names(x)==sampleid_fecal);return(x[idx])})))
    Y=rbind(Y,as.numeric(sapply(oraldata,function(x){idx=(names(x)==sampleid_oral);return(x[idx])})))
  }

  #remove NA in X and Y
  xna=which(apply(X,1,function(x){is.na(sum(x))})==T)
  yna=which(apply(Y,1,function(x){is.na(sum(x))})==T)
  xyna=union(xna,yna)
  X=X[-xyna,]
  Y=Y[-xyna,]


  #obtain the taxa
  OTUname_fecal = sapply(AG_Taxa[OTU_fecal],function(x)x$id)
  colnames(X)=OTUname_fecal
  Taxa_fecal = t(as.data.frame( lapply(AG_Taxa[OTU_fecal],function(x) x$metadata$taxonomy))); rownames(Taxa_fecal) = OTUname_fecal
  for(i in 1:nrow(Taxa_fecal))
    for(j in 1:ncol(Taxa_fecal)) Taxa_fecal[i,j] = substring(Taxa_fecal[i,j], 4,1000)

  OTUname_oral = sapply(AG_Taxa[OTU_oral],function(x)x$id)
  colnames(Y)=OTUname_oral
  Taxa_oral = t(as.data.frame( lapply(AG_Taxa[OTU_oral],function(x) x$metadata$taxonomy))); rownames(Taxa_oral) = OTUname_oral
  for(i in 1:nrow(Taxa_oral))
    for(j in 1:ncol(Taxa_oral)) Taxa_oral[i,j] = substring(Taxa_oral[i,j], 4,1000)

  #obtain the tree
  Tree_fecal=prune_taxa(AG_tree$tip.label %in%  OTUname_fecal, AG_tree)
  Tree_oral=prune_taxa(AG_tree$tip.label %in%  OTUname_oral, AG_tree)
  #sort X column s.t. it matches the order of the tree
  X=X[,order(match(colnames(X),Tree_fecal$tip.label))]
  Y=Y[,order(match(colnames(Y),Tree_oral$tip.label))]

  return(list(X=X,Y=Y,Tree_S1=Tree_fecal,Tree_S2=Tree_oral,Taxa_S1=Taxa_fecal,Taxa_S2=Taxa_oral))
}
