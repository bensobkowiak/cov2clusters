numberClusters<-function(acceptTrans){
  names<-unique(c(as.character(acceptTrans[,1]),as.character(acceptTrans[,2])))
  cluster_results<-matrix(NA,ncol = 2, nrow = length(names))
  cluster_results[,1]<-names
  if (nrow(cluster_results)>0){
    nums<-1:nrow(cluster_results)
    for (i in 1:nrow(cluster_results)){
      if (i%in%nums){
        clustnums<-which(cluster_results[,1] %in% unique(c(acceptTrans[which(acceptTrans[,1]==cluster_results[i,1] | 
                                                                               acceptTrans[,2]==cluster_results[i,1]),1:2])))
        prevClust<-unique(cluster_results[clustnums,2])
        prevClust<-as.numeric(prevClust[!is.na(prevClust)])
        if (length(prevClust)>0){
          cluster_results[unique(c(which(cluster_results[,2]%in%prevClust),clustnums)),2]<-min(prevClust)
        } else{
          cluster_results[clustnums,2]<-min(nums)
        }
        nums<-nums[!nums%in%clustnums]
      }
    }
  }
  return(cluster_results)
}
