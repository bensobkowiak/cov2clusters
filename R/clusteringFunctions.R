
### Logit probability function
TransProbs<-function(dataInput,contactData,probThreshold,dateThreshold,restrictClusters,beta){
  PatDist<-dataInput$PatDist
  dates<-dataInput$dates
  if (contactData){
    contacts<-dataInput$contacts
  }
  if (restrictClusters){
    variable<-dataInput$variable
  }
  
  transmat<-foreach(i=seq(ncol(PatDist)-1), .combine = rbind) %dopar% {
    transmatrix<-matrix(ncol = 3)
    for (j in seq(i+1, ncol(PatDist))){
      if (contactData){
        covariates<-c(1,0,0,0) #Empty covariates
        covariates[2]<-PatDist[i,j]
        covariates[3]<-abs(dates[i]-dates[j])
        covariates[4]<-contacts[i,j]
      } else {
        covariates<-c(1,0,0) #Empty covariates
        covariates[2]<-PatDist[i,j]
        covariates[3]<-abs(dates[i]-dates[j])
        beta<-beta[1:3]
      }
      
      # Run logit if date threshold met
      if (!restrictClusters | restrictClusters && variable[i]==variable[j]){
        if (!is.na(dateThreshold) && covariates[3]>(dateThreshold)){
          trans<-0
        } else {
          trans = 1/(1 + exp(- (sum(beta*covariates))))
        }
      }
      # Append to transmatrix if > min probThreshold
      if (trans>=min(probThreshold)){
        transres<-c(colnames(PatDist)[i],colnames(PatDist)[j],trans)
        transmatrix<-rbind(transmatrix,transres)
      }
    }
    transmatrix
  }
  colnames(transmat)<-c("host1","host2","Prob")
  transmat<-transmat[!is.na(transmat[,3]),]
}

### Number cluster function
numberClusters<-function(acceptTrans){
  names<-unique(c(as.character(acceptTrans[,1]),as.character(acceptTrans[,2])))
  cluster_results<-matrix(NA,ncol = 2, nrow = length(names))
  cluster_results[,1]<-names
  if (nrow(cluster_results)>0){
    nums<-1:nrow(cluster_results)
    for (i in 1:nrow(cluster_results)){
      if (i%in%nums){
        clustnums<-which(cluster_results[,1] %in% unique(unlist(c(acceptTrans[which(acceptTrans[,1]==cluster_results[i,1] | 
                                                                                      acceptTrans[,2]==cluster_results[i,1]),1:2]))))
        assocnums<-which(cluster_results[,1] %in% unique(unlist(c(acceptTrans[which(acceptTrans[,1] %in% cluster_results[clustnums,1] | 
                                                                                      acceptTrans[,2] %in% cluster_results[clustnums,1]),1:2]))))
        allnums<-c(clustnums,assocnums)
        prevClust<-unique(cluster_results[allnums,2])
        prevClust<-as.numeric(prevClust[!is.na(prevClust)])
        if (length(prevClust)>0){
          cluster_results[unique(c(which(cluster_results[,2]%in%prevClust),allnums)),2]<-min(prevClust)
        } else{
          cluster_results[allnums,2]<-min(nums)
        }
        nums<-nums[!nums%in%clustnums]
      }
    }
  }
  return(cluster_results)
}
