#############################################################################################
########### cov2clusters - SARS-CoV-2 genomic clusters from phylogenetic trees ##############
#############################################################################################

#' @param treeName Name of the input phylogenetic tree, in newick format
#' @param metafile Name of the date file if in .csv format with columns "sequence ID" and "dates" (not required if json_dates = TRUE and restrictClusters = FALSE)
#' @param json_dates If TRUE - dates are supplied in json format file
#' @param json_file If json_dates = TRUE - name .json dates file
#' @param contactData If TRUE - file will be supplied with contact data 
#' @param contactFile File contacting contact, in .csv format (required if contactData = TRUE)
#' @param beta List of the beta coefficients beta0 - beta3 (Intercept, patristic distance, date difference, and contact data(if applicable))
#' @param returnTransProbs Return text file of pairwise logit probabilities (must = TRUE if newClustering = FALSE)
#' @param dateThreshold Integer for the hard upper limit of date difference (in days) to link pairs of sequences in a cluster (default = 40)
#' @param restrictClusters If TRUE - Cluster sequences only with a shared variable in the metafile column 'restrictCluster' (optional)
#' @param probThreshold Integer or list - Pairwise probability threshold to cluster sequences
#' @param newClustering If TRUE - 
#' @param pastTransProbs If TRUE - 
#' @param clusterFile If TRUE - 
#' @param clusternameIdent If TRUE - 
#' @param outfile Prefix for the output file
#' @param no.Cores Integer - number of cores to run the function, can be multithreaded
#' @return Text file with three columns - Sequence ID from tree tip, cluster name, and past cluster name (if applicable)
#' @export

cov2clusters<-function(treeName="tree.nwk",metafile=NA,
                       json_dates=TRUE,json_file="branch_lengths.json",
                       contactData=FALSE,contactFile=NA,
                       beta=c(3,-19735.98,-0.075,-0.2),returnTransProbs=FALSE,
                       dateThreshold=40,restrictClusters=FALSE,
                       probThreshold=c(0.8,0.9),
                       newClustering=TRUE, pastTransProbs=NA,
                       clusterFile=NA,
                       clusternameIdent = "clust",
                       outfile="SARS-CoV-2",no.Cores=1){
  # Load packages
  require(ape)
  require(reshape2)
  require(rjson)
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
  dataInput$PatDist<-dataInput$PatDist[order(row.names(dataInput$PatDist)),order(colnames(dataInput$PatDist))]
  
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
  # remove patristic distance without dates
  dataInput$PatDist<-dataInput$PatDist[which(colnames(dataInput$PatDist) %in% dates[,1]),
                                       which(colnames(dataInput$PatDist) %in% dates[,1])]
  dataInput$dates<-round(lubridate::decimal_date(as.Date(dates[,2]))*365)
  
  if (restrictClusters){
    variable <- read.csv(metafile,check.names = F)
    variable<-variable[which(variable[,1] %in% colnames(dataInput$PatDist)),]
    variable<-variable[order(variable[,1]),]
    dataInput$variable<-variable$restrictCluster
  }
  # Contact data
  if (contactData){
    contacts<-read.table(contactFile)
    contacts<-contacts[which(row.names(contacts) %in% colnames(dataInput$PatDist)),
                       which(colnames(contacts) %in% colnames(dataInput$PatDist))]
    contacts<-contacts[order(colnames(contacts)),order(colnames(contacts))]
    dataInput$contacts<-contacts
  }
  
  # Transmission matrix
  transmat<-TransProbs(dataInput,contactData,probThreshold,dateThreshold,restrictClusters,beta)
   if (!newClustering){
    pastTransProbs<-read.table(pastTransProbs,header = T,check.names = F)
    newnames<-colnames(dataInput$PatDist)[-which(colnames(dataInput$PatDist) %in% 
                                                  pastTransProbs$host1 |
                                                   colnames(dataInput$PatDist) %in% 
                                                   pastTransProbs$host2)]
    transmat<-transmat[which(transmat[,1] %in% newnames | 
                               transmat[,2] %in% newnames),]
    if (nrow(transmat)>0){
      transmat<-rbind(transmat,pastTransProbs)
    } else {
      transmat<-pastTransProbs
    }
  }
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
      clustersDF<-data.frame(ClusterName=rep(NA,length(clusters)),ClusterSize=NA,ClusterComposition=rep(NA,length(clusters)))
      months<-c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
      randstring<-stri_rand_strings(length(clusters)*3,3,"[b-df-hj-np-tv-z]")
      randstring<-unique(randstring) # make sure string is unique
      pastClusterString<-stri_sub(unique(pastClusters$Cluster_no[!pastClusters$Cluster_no=="-1"]),-3)
      randstring<-randstring[!randstring%in%pastClusterString] # not in past cluster string
      for (clust in 1:length(clusters)){
        rownum<-which(cluster_results[,2]==clusters[clust])
        clustersDF$ClusterSize[clust]<-length(rownum)
        pastClusterNames<-pastClusters[which(pastClusters[,1] %in% cluster_results[rownum,1]),2]
        for (j in 1:length(rownum)){
          if (cluster_results[rownum[j],1] %in% pastClusters[,1]){
            cluster_results[rownum[j],3]<-pastClusters[which(pastClusters[,1] %in% cluster_results[rownum[j],1]),2]
          }}
        if (length(pastClusterNames)>0){
          pastClusterDF<-data.frame(table(pastClusterNames),prop.in.clust=0)
          for (i in 1:nrow(pastClusterDF)){
            if (pastClusterDF$pastClusterNames[i]=="-1"){
              pastClusterDF$prop.in.clust[i]<-NA
            } else {
              pastClusterDF$prop.in.clust[i]<-round(pastClusterDF$Freq[i]/length(which(pastClusters[,2] %in% pastClusterDF$pastClusterNames[i])),2)
            }}
          clustersDF$ClusterComposition[clust]<-paste0(paste0(pastClusterDF$pastClusterNames,",",pastClusterDF$Freq,",",pastClusterDF$prop.in.clust),":",collapse = "")
          pastClusterDF<-pastClusterDF[pastClusterDF$pastClusterNames!="-1",]
          if (nrow(pastClusterDF)>0 &&
              pastClusterDF[which.max(pastClusterDF$Freq),3]>=0.6){
            clustername<-as.character(pastClusterDF[which.max(pastClusterDF$Freq),1])
            cluster_results[rownum,2]<-clustername
            clustersDF$ClusterName[clust]<-clustername
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
    colnames(cluster_results)<-c("SampleID","Cluster.No","pastCluster.No")
    if (!newClustering){
      pastnoclust<-cbind(as.character(pastClusters[which(!as.character(pastClusters$SampleID) %in% 
                                                           cluster_results[,1]),1]),"-1",
                         as.character(pastClusters[which(!as.character(pastClusters$SampleID) %in% 
                                                           cluster_results[,1]),2]))
      if (ncol(pastnoclust)==3){
        cluster_results<-rbind(as.matrix(cluster_results),pastnoclust)
      }
    }
    noclust<-cbind(as.character(dates[which(!as.character(dates[,1]) %in% cluster_results[,1]),1]),"-1",NA)
    if (ncol(noclust)==3){
      cluster_results<-rbind(as.matrix(cluster_results),noclust)
    }
    write.table(cluster_results,paste0(outfile,"_",probThreshold[threshold],"_GenomicClusters",Sys.Date(),".txt"),row.names = F,sep = "\t",quote = F)
  }
}


