## Main clustering function

cov2clusters<-function(treeName="tree.nwk",metafile=NA,
                       json_dates= TRUE,json_file="branch_lengths.json",
                       beta=c(3,-19735.98,-0.075),returnTransProbs=FALSE,
                       dateThreshold=40,restrictClusters=FALSE,
                       probThreshold=c(0.7,0.8,0.9),
                       newClustering=TRUE, clusterFile=NA,
                       clusternameIdent = "clust",
                       outfile="SARS-CoV-2",no.Cores=1){
  # Load packages
  require(ape)
  require(reshape2)
  require(rjson)
  require(dplyr)
  require(lubridate)
  require(stringi)
  require(foreach)
  require(doMC)
  registerDoMC(no.Cores)
  options(stringsAsFactors = F)
  
  #Data Setup
  dataInput<-list()
  
  #Patristic distance
  tree<-read.tree(treeName)
  dataInput$PatDist<-as.matrix(cophenetic.phylo(tree))
  dataInput$PatDist<-dataInput$PatDist[order(colnames(dataInput$PatDist)),order(row.names(dataInput$PatDist))]
  
  #Load Metadata
  if (json_dates){
    json_data <- fromJSON(file=json_file)
    dates<-data.frame(names=colnames(dataInput$PatDist),date=NA)
    for (i in 1:nrow(dates)){
      dates$date[i]<-as.character(unlist(json_data$nodes[[which(names(json_data$nodes)==dates$names[i])]][3]))
    }
  } else {
    dates <- read.csv(metafile,check.names = F)
  }
  dates<-dates[which(dates[,1] %in% colnames(dataInput$PatDist)),]
  dates<-dates[order(dates[,1]),]
  dataInput$dates<-lubridate::decimal_date(as.Date(dates[,2]))*365
  
  if (restrictClusters){
    variable <- read.csv(metafile,check.names = F)
    variable<-variable[which(variable[,1] %in% colnames(dataInput$PatDist)),]
    variable<-variable[order(variable[,1]),]
    dataInput$variable<-variable[,3]
  }
  transmat<-TransProbs(dataInput,probThreshold,dateThreshold,restrictClusters,beta)
  if (returnTransProbs){
    write.table(transmat,paste0(outfile,"_TransProbs_",Sys.Date(),".txt"),quote = F,sep = "\t",row.names = F)
  }
  
  # Cluster using defined thresholds 
  for (threshold in 1:length(probThreshold)){
    acceptTrans<-transmat[which(transmat[,3]>=probThreshold[threshold]),]
    
    # Number clusters
    cluster_results<-numberClusters(acceptTrans)
    cluster_results<-data.frame(cluster_results,pastClustering=NA)
    ## Rename clusters by past cluster names or by earliest case date and location if new clustering
    if (newClustering){
      clusters<-unique(cluster_results[,2])
      months<-c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
      randstring<-stri_rand_strings(length(clusters)*2,3,"[b-df-hj-np-tv-z]")
      randstring<-unique(randstring) # make sure string is unique
      for (clust in 1:length(clusters)){
        mon<-month(min(dates[which(dates[,1] %in% 
                                     cluster_results[which(cluster_results[,2]==clusters[clust]),1]),2]))
        yea<-substr(year(min(dates[which(dates[,1] %in% 
                                           cluster_results[which(cluster_results[,2]==clusters[clust]),1]),2])),3,4)
        cluster_results[which(cluster_results[,2]==clusters[clust]),2]<-paste0(paste0(months[mon],yea),".",clusternameIdent,".",randstring[clust])
      }
    } else {
      pastClusters<-read.table(clusterFile,header = T,check.names = F)
      clusters<-unique(cluster_results[,2])
      clustersDF<-data.frame(ClusterName=rep(NA,length(clusters)),ClusterComposition=rep(NA,length(clusters)))
      months<-c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
      randstring<-stri_rand_strings(length(clusters)*2,3,"[b-df-hj-np-tv-z]")
      randstring<-unique(randstring) # make sure string is unique
      for (clust in 1:length(clusters)){
        rownum<-which(cluster_results[,2]==clusters[clust])
        pastClusterNames<-pastClusters[which(pastClusters[,1] %in% cluster_results[rownum,1]),2]
        if (length(pastClusterNames)>0){
          pastClusterDF<-data.frame(table(pastClusterNames),prop.in.clust=0)
          for (i in 1:nrow(pastClusterDF)){
            pastClusterDF$prop.in.clust[i]<-pastClusterDF$Freq/length(which(pastClusters[,2] %in% pastClusterDF$pastClusterNames[i]))
          }
          clustersDF$ClusterComposition[clust]<-paste0(paste0(pastClusterDF$pastClusterNames,",",pastClusterDF$prop.in.clust),",",collapse = "")
          pastClusterDF<-pastClusterDF[pastClusterDF$pastClusterNames!="-1",]
          if (nrow(pastClusterDF)>0 &
              pastClusterDF[which.max(pastClusterDF$Freq),3]>=0.6){
            clustername<-as.character(pastClusterDF[which.max(pastClusterDF$Freq),1])
            cluster_results[rownum,2]<-clustername
            clustersDF$ClusterName[clust]<-clustername
            cluster_results[rownum,3]<-pastClusterNames
          } else {
            mon<-month(min(dates[which(dates[,1] %in% 
                                         cluster_results[rownum,1]),2]))
            yea<-substr(year(min(dates[which(dates[,1] %in% 
                                               cluster_results[rownum,1]),2])),3,4)
            clustername<-paste0(paste0(months[mon],yea),".",clusternameIdent,".",randstring[clust])
            cluster_results[rownum,2]<-clustername
            clustersDF$ClusterName[clust]<-clustername
          }
        } else {
          mon<-month(min(dates[which(dates[,1] %in% 
                                       cluster_results[rownum,1]),2]))
          yea<-substr(year(min(dates[which(dates[,1] %in% 
                                             cluster_results[rownum,1]),2])),3,4)
          clustername<-paste0(paste0(months[mon],yea),".",clusternameIdent,".",randstring[clust])
          cluster_results[rownum,2]<-clustername
          clustersDF$ClusterName[clust]<-clustername
        }
      }
      write.csv(clustersDF,paste0(outfile,"_",probThreshold[threshold],"_ClustersSummary",Sys.Date(),".csv"),row.names = F)# cluster summary output
    }
    
    # Output file 
    cluster_results<-cluster_results[order(cluster_results[,2]),]
    row.names(cluster_results)<-NULL
    colnames(cluster_results)<-c("SampleID","Cluster_no","pastClustering")
    noclust<-cbind(as.character(dates[which(!as.character(dates[,1]) %in% cluster_results[,1]),1]),"-1",NA)
    if (!newClustering){
      for (i in 1:nrow(noclust)){
        if (is.element(noclust[i,1],pastClusters[,1])){
          noclust[i,3]<-pastClusters[which(pastClusters[,1]==noclust[i,1]),2]
        }
      }
    }
    if (ncol(noclust)>1){
      cluster_results<-rbind(as.matrix(cluster_results),noclust)
    }
    write.table(cluster_results,paste0(outfile,"_",probThreshold[threshold],"_GenomicClusters",Sys.Date(),".txt"),row.names = F,sep = "\t",quote = F)
  }
}
