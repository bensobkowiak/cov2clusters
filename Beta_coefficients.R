
covidBeta<-function(metafile,distanceMatrix){
  
  
  metaData <- read.csv(metafile)
  distMatrix<-read.table(distanceMatrix)
  
  #order metaData by distMatrix
  metaData<-metaData[which(metaData[,1] %in% colnames(distMatrix)),]
  meta.order<-integer()
  for (i in 1:nrow(metaData)){
    meta.order[i]<-which(colnames(distMatrix)==metaData[i,1])
  }
  metaData<-metaData[order(meta.order),]
  
  ## Remove no epi
  newdistMatrix<-distMatrix[-which(is.na(metaData$Cluster_ID)),-which(is.na(metaData$Cluster_ID))]
  newmetaData<-metaData[-which(is.na(metaData$Cluster_ID)),]
  
  # Set up input data list
  dataInput<-list()
  dataInput$IDs<-colnames(newdistMatrix)
  dataInput$dist<-newdistMatrix
  dataInput$dates<-lubridate::decimal_date(as.Date(newmetaData[,2]))
  transres<-matrix(NA,ncol = ncol(newdistMatrix),nrow = nrow(newdistMatrix))
  for (i in seq(length(dataInput$IDs)-1)){
    for (j in seq(i+1, length(dataInput$IDs))){
      if (newmetaData$Cluster_ID[which(newmetaData$SamID==dataInput$IDs[i])] == newmetaData$Cluster_ID[which(newmetaData$SamID==dataInput$IDs[j])]){
        transres[i,j] <- 1
      } else {
        transres[i,j] <- 0
      }
    }
  }
  colnames(transres)<-dataInput$IDs
  row.names(transres)<-dataInput$IDs
  transmat<-reshape2::melt(transres, varnames = c('row', 'col'), na.rm = TRUE)
  
  dateres<-matrix(NA,ncol = ncol(newdistMatrix),nrow = nrow(newdistMatrix))
  for (i in seq(length(dataInput$IDs)-1)){
    for (j in seq(i+1, length(dataInput$IDs))){
      dateres[i,j] <- abs(dataInput$dates[i]-dataInput$dates[j])
    }
  }
  colnames(dateres)<-dataInput$IDs
  row.names(dateres)<-dataInput$IDs
  datemat<-reshape2::melt(dateres, varnames = c('row', 'col'), na.rm = TRUE)
  row.names(transmat)<-NULL
  row.names(datemat)<-NULL
  newdistMatrix[lower.tri(newdistMatrix,diag = T)] = NA
  distmat<-reshape2::melt(as.matrix(newdistMatrix), varnames = c('row', 'col'), na.rm = TRUE)
  row.names(distmat)<-NULL
  
  # Paste together 
  Mat4reg<-cbind(distmat,datemat$value,transmat$value)
  Mat4reg<-Mat4reg[-which(Mat4reg$`datemat$value`>(21/365)),]
  
  table(Mat4reg$trans)
  colnames(Mat4reg)<-c("sample1","sample2","distance","date","trans")
  
  # do log regression
  res<-glm(trans ~ log(distance+1) + date,data = Mat4reg,family = 'binomial') 
  return(as.numeric(summary(res)$coef[,1]))
}
