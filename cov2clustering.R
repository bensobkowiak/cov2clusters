
cov2clusters <- function(metafile, distanceMatrix, beta = c(1, 1, 1),
                         SNPthreshold = NULL, includeDate = TRUE, dateThreshold=NULL,  
                         includeSpatial = FALSE, CSVfile=FALSE,probThreshold=0.8){
  
  #Check inputs
  if (((includeDate & !includeSpatial) |(!includeDate & includeSpatial))  && length(beta)!=3){
    print("Error: not enough beta coefficients supplied")
    return()
  } else if ((includeDate && includeSpatial) && length(beta)!=4){
    print("Error: not enough beta coefficients supplied")
    return()
  }
  
  #Load data
  metaData <- read.csv(metafile)
  distMatrix<-read.table(distanceMatrix)

  #order metaData by distMatrix
  metaData<-metaData[which(metaData[,1] %in% colnames(distMatrix)),]
  meta.order<-integer()
  for (i in 1:nrow(metaData)){
    meta.order[i]<-which(colnames(distMatrix)==metaData[i,1])
  }
  metaData<-metaData[order(meta.order),]
  
  # Set up input data list
  dataInput<-list()
  dataInput$IDs<-colnames(distMatrix)
  dataInput$dist<-distMatrix
  if (includeDate){dataInput$dates<-lubridate::decimal_date(as.Date(metaData[,2]))}
  if (includeSpatial){dataInput$spatial<-metaData[,3]}
  dataInput$trans<-matrix(NA,ncol = ncol(distMatrix),nrow = nrow(distMatrix))
  
  # Calculate prob of transmission
  for (i in seq(length(dataInput$IDs)-1)){
    for (j in seq(i+1, length(dataInput$IDs))){
      covariates<-integer()
      covariates[1]<-1
      covariates[2]<-log(dataInput$dist[i,j]+1)
      if (includeDate && !includeSpatial){
        covariates[3]<-abs(dataInput$dates[i]-dataInput$dates[j])
      } else if (!includeDate && includeSpatial){
        if (dataInput$spatial[i]==dataInput$spatial[j]){
          covariates[3]<-1
        } else {
          covariates[3]<-0
        }
      } else if (includeDate && includeSpatial){
        covariates[3]<-abs(dataInput$dates[i]-dataInput$dates[j])
        if (dataInput$spatial[i]==dataInput$spatial[j]){
          covariates[4]<-1
        } else {
          covariates[4]<-0
        }
      }
      if (!is.null(SNPthreshold) && covariates[2]>=SNPthreshold){
        dataInput$trans[i,j]<-0
      } else if (!is.null(dateThreshold) && covariates[3]>(dateThreshold/365)){
        dataInput$trans[i,j]<-0
      #} else if (covariates[2]==0){
      #  dataInput$trans[i,j]<-1
      } else {
        dataInput$trans[i,j] = 1/(1 + exp(- (sum(beta*covariates))))
      }
    }
  }
  colnames(dataInput$trans)<-dataInput$IDs
  row.names(dataInput$trans)<-dataInput$IDs
  
  # Accept probs >= probThreshold
  transmat<-reshape2::melt(dataInput$trans, varnames = c('row', 'col'), na.rm = TRUE)
  acceptTrans<-transmat[which(transmat$value>=probThreshold),]
  
  # Assign clusters
  names<-unique(c(as.character(acceptTrans[,1]),as.character(acceptTrans[,2])))
  cluster_results<-as.data.frame(matrix(0,ncol = 2, nrow = length(names)))
  cluster_results[,1]<-names
  nums<-1:nrow(acceptTrans)
  colnames(acceptTrans)<-c("sample1","sample2","value")
  if (nrow(cluster_results)>0){
    for (i in 1:nrow(acceptTrans)){
      if (cluster_results[which(cluster_results[,1]==acceptTrans$sample1[i]),2]==0 &cluster_results[which(cluster_results[,1]==acceptTrans$sample2[i]),2]==0){
        cluster_results[which(cluster_results[,1]==acceptTrans$sample1[i]),2]<-min(nums)
        cluster_results[which(cluster_results[,1]==acceptTrans$sample2[i]),2]<-min(nums)
        nums<-nums[-1]
      } else {
        a<-cluster_results[which(cluster_results[,1]==acceptTrans$sample1[i]),2]
        b<-cluster_results[which(cluster_results[,1]==acceptTrans$sample2[i]),2]
        if (a!=b){
          c<-c(a,b)
          e<-which(c==0)
          if (length(e)>0){
            c<-c[-e]
          }
          if (length(c)==1){
            cluster_results[which(cluster_results[,1]==acceptTrans$sample1[i]),2]<-c
            cluster_results[which(cluster_results[,1]==acceptTrans$sample2[i]),2]<-c
          } else {
            cluster_results[which(cluster_results[,1]==acceptTrans$sample1[i]),2]<-min(c)
            cluster_results[which(cluster_results[,1]==acceptTrans$sample2[i]),2]<-min(c)
            d<-which(cluster_results[,2]==max(c))
            cluster_results[d,2]<-min(c)
          }
        }
      }
    }
  }
  
  ## Rename clusters by earliest case date and location
  
  
  # Output file
  cluster_results<-cluster_results[order(cluster_results[,2]),]
  row.names(cluster_results)<-NULL
  colnames(cluster_results)<-c("ID","Cluster_no")
  noclust<-cbind(as.character(metaData[which(!as.character(metaData[,1]) %in% cluster_results$ID),1]),NA)
  cluster_results<-rbind(as.matrix(cluster_results),noclust)
  if (CSVfile){
  write.csv(cluster_results,"cov2cluster_results.csv",row.names = F)
  }
  return(cluster_results)
} 
