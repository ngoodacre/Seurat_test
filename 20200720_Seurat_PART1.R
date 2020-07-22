####################################################################################################
#### The following is a rendition of Shane Lofgren's Seurat/Monocle pipeline
#### based on a 4-page written summary and using functions mainly from the tutorial:
#### https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
#### as well as other tutorials from:
#### https://satijalab.org/seurat/vignettes.html
####################################################################################################


######################################>->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
######################################      RUN PART 1      
######################################>->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->

projectdir<-getwd()
scriptsdir<-paste(projectdir,"scripts",sep="/")
source(paste(scriptsdir,"20200721_Seurat_fxns_PART0.R",sep="/"))#load supporting fxns
source(paste(scriptsdir,"20200720_Seurat_fxns_PART1.R",sep="/"))#load supporting fxns
source(paste(scriptsdir,"20200721_Seurat_fxns_GENERIC.R",sep="/"))#load supporting fxns
datadir<-paste(projectdir,"h5",sep="/")
refdir<-paste(projectdir,"ref",sep="/")
args = commandArgs(trailingOnly=TRUE)
project=args[1]
date=args[2]
analysis=args[3]
outdir=paste(projectdir,"Out",sep="/")
sessionsdir=paste(projectdir,"Sessions",sep="/")
figs_folder=paste(outdir,"Figs",sep="/")
ngroup_levels=as.integer(args[4])
splitby=args[5]
refgenes_file<-paste(refdir,args[6],sep="/")
refgenes<-unique(read.delim(refgenes_file,sep='\n')[[1]])
setup_project_env(projectdir)
badfiles=strsplit(args[7],',')


session_name1<-paste(sessionsdir,paste(project,date,analysis,"Seurat_PART1.rds",sep="_"),sep="/")
message(paste("Session will be saved as ",session_name1,sep=""))
scrna.bysplit<-RUN_PART_1(datadir,badfiles,refgenes,ngroup_levels,splitby,session_name1,figs_folder,date)
#####->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
##