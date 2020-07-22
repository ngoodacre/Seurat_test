####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
#######################################       PROJECT SET-UP        ################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(viridis)
library(pals)
setup_project_env<-function(project_folder){
  print(paste("Initializing project ",project_folder,sep=""))
  out_folder=paste(project_folder,"Out",sep="/")
  sessions_folder=paste(project_folder,"Sessions",sep="/")
  figs_folder=paste(out_folder,"Figs",sep="/")
  if(dir.exists(project_folder)!=TRUE){
    dir.create(project_folder)
    print(paste("Creating ",project_folder,sep=""))
  }
  if(dir.exists(out_folder)!=TRUE){
    dir.create(out_folder)
    print(paste("Creating ",out_folder,sep=""))
  }
  if(dir.exists(sessions_folder)!=TRUE){
    dir.create(sessions_folder)
    print(paste("Creating ",sessions_folder,sep=""))
  }
  if(dir.exists(figs_folder)!=TRUE){
    dir.create(figs_folder)
    print(paste("Creating ",figs_folder,sep=""))
  }
}