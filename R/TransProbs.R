# Logit probability function
TransProbs<-function(dataInput,probThreshold,dateThreshold,restrictClusters,beta){
  PatDist<-dataInput$PatDist
  dates<-dataInput$dates
  if (restrictClusters){
    variable<-dataInput$variable
  }
  
  transmat<-foreach(i=seq(ncol(PatDist)-1), .combine = rbind) %dopar% {
    transmatrix<-matrix(ncol = 3)
    for (j in seq(i+1, ncol(PatDist))){
      
      covariates<-c(1,0,0) #Empty covariates
      covariates[2]<-PatDist[i,j]
      covariates[3]<-abs(dates[i]-dates[j])
      
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
