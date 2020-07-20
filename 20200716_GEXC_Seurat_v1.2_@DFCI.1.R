####################################################################################################
#### The following is a rendition of Shane Lofgren's Seurat/Monocle pipeline
#### based on a 4-page written summary and using functions mianly from the tutorial:
#### https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
#### as well as other tutorials from:
#### https://satijalab.org/seurat/vignettes.html
####################################################################################################





##
#####->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(viridis)
source("/Users/ngoodacre/Desktop/pipeline_protocols/GEX_Seurat/pipeline/v1/20200707_Seurat_fxns.R")
#####->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
##





####################################################################################################
#### PART 0 - INITIALIZE PROJECT, PATHING


#####->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
setwd("/Users/ngoodacre/Desktop/pipeline_protocols/GEX_Seurat/3Tdata/DFCI.1/h5")
ref_folder='/Users/ngoodacre/Desktop/pipeline_protocols/GEX_Seurat/3Tdata/DFCI.1/ref'
project_folder='/Users/ngoodacre/Desktop/pipeline_protocols/GEX_Seurat/3Tanalyses/DFCI.1_20200715'
out_folder=paste(project_folder,"Out",sep="/")
sessions_folder=paste(project_folder,"Sessions",sep="/")
figs_folder=paste(out_folder,"Figs",sep="/")
analysis='DFCI.1'
date='20200715-p'
#####->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
##


##
#### Only run below if first time
r<-setup_project_env(project_folder,getwd())
num_files<-r[[1]]
ngroup_levels<-r[[2]]
####
##


####################################################################################################










#########################################################################################################################
#########################################################################################################################
######################################      PART I -> FINAL GEX MATRIX      #############################################
#########################################################################################################################
#########################################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

RUN_PART_1<-function(bad_files,nfeatures,refgenes,npcas,ngroup_levels,session_name){
  h5files=list.files(getwd())
  r<-setup_seurat_objects(getwd(),h5files,bad_files,ngroup_levels,analysis)
  seurat_objects<-r[[1]]
  samples<-r[[2]]
  seurat_objects<-patch_gene_rows(seurat_objects,refgenes,analysis)
  scrna.merged<-merge_seurat_objects_DFCI.1(seurat_objects,samples,analysis,ngroup_levels)
  scrna.bypatient<-split_merged_seurat_object_bypatient(scrna.merged,splitby,figs_folder,date)
  scrna.integrated<-integrate_seurat_bypatient(scrna.bypatient,nfeatures,refgenes,npcas,ngroup_levels)
  #### Save x1 
  saveRDS(scrna.integrated, file = paste(sessions_folder,session_name,sep="/"))
  return(scrna.integrated)
}


##
#####->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
bad_files=c()
# DEFAULT PARAMS BELOW
nfeatures=2000
npcas=30
#
refgenes_file<-paste(ref_folder,'20200714_combined_markers.txt',sep="/")
refgenes<-unique(read.delim(ref_genes_file,sep='\n')[[1]])
ngroup_levels=3
splitby='patient_timepoint'
#splitby='patient'
date="20200715-pt"
session_name<-paste(analysis,"_",date,"_","Seurat_pt1.rds",sep="")
scrna.integrated<-RUN_PART_1(bad_files,nfeatures,refgenes,npcas,ngroup_levels,session_name)
#####->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
##





#########################################################################################################################
#########################################################################################################################
######################################      PART II -> DIM REDUC, UMAP, CLUSTERING      #################################
#########################################################################################################################
#########################################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

RUN_PART_2<-function(session_name1,sessions_name2,npcas,res,figs_folder_global,groups,date){
  ####################################################################################################
  #### LOAD THE  DATA FROM PART 1
  scrna.integrated=readRDS(paste(sessions_folder,session_name1,sep="/"))
  ####################################################################################################
  #### SET UP OUTPUT FOLDER STRUCTURE
  setup_figs_outs(scrna.integrated,figs_folder)
  ####################################################################################################
  #### SCALING THE DATA
  all.genes <- rownames(scrna.integrated)
  scrna.integrated <- ScaleData(scrna.integrated, features = all.genes)
  ####################################################################################################
  #### PERFORMING DIMENSIONALITY  EXPLORATION AND REDUCTION
  # DEFAULT PARAM
  # npcas=30
  scrna.integrated <- explore_dimensionality(scrna.integrated,npcas,figs_folder_global,date)
  ####################################################################################################
  #### CLUSTERING CELLS (LOUVAIN ALGORITHM)
  # DEFAULT PARAM
  # res=0.5
  scrna.clustered<-cluster_cells(scrna.integrated,npcas,res)
  ####################################################################################################
  #### SET UP CLUSTERED OBJECT FIGS OUTS FOLDER
  setup_figs_outs_clustered(scrna.clustered,figs_folder)
  ####################################################################################################
  #### PLOTTING GLOBAL UMAPS (by patient, by sort, by sample (patient+sort))
  plot_global_umaps(scrna.clustered,figs_folder_global,groups,date)
  #### Save x2 
  saveRDS(scrna.clustered, file = paste(sessions_folder,session_name2,sep="/"))
  return(scrna.clustered)
}


##
#####->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
npcas<-30
res<-0.5
session_name1<-paste(analysis,"_",date,"_","Seurat_pt1.rds",sep="")
session_name2<-paste(analysis,"_",date,"_","Seurat_pt2.rds",sep="")
figs_folder_global<-paste(figs_folder,"global",sep="/")
groups=c('orig.ident','sort','patient','timepoint','patient_timepoint','seurat_clusters')
scrna.clustered<-RUN_PART_2(session_name1,session_name2,npcas,res,figs_folder_global,date)
#####->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
##





#########################################################################################################################
#########################################################################################################################
###################################      PART III -> DIFFERENTIAL EXPRESSION ANALYSIS     ###############################
#########################################################################################################################
#########################################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
RUN_PART_3<-function(session_name2,session_name3,figout,dimreduce_method,minpct,logfcthreshold,bm_file1,bm_file2){
  ####################################################################################################
  #### LOAD THE DATA FROM PART 2
  scrna.clustered=readRDS(paste(sessions_folder,session_name2,sep="/"))
  signs<-parse_celltype_signatures(celltype_signature_ref)
  #setup_figs_outs_celltypes(signs,figs_folder)
  r<-define_cell_types(scrna.clustered,signs)
  scrna.clustered<-r[[1]]
  cell_types<-r[[2]]
  
  plot_cell_types__main(scrna.clustered,signs,figs_folder,date)
  
  
    
  #visualize_biomarkers(scrna.clustered,figs_folder,dimreduce_method,bm_file1,bm_file2,date)
  
  visualize_biomarkers_dotplot(scrna.clustered,figs_folder,dimreduce_method,bm_file1,bm_file2,date)
  r<-find_clustermarkers(scrna.clustered,figout,dimreduce_method,minpct,logfcthreshold,figs_folder_clustered,date)#line ~367
    
  scrna.de<-r[0]
  cluster_markers_list<-r[1]
  #### Save x3 
  saveRDS(scrna.de, file = paste(sessions_folder,session_name3,sep="/"))
  return(r)
}


##
#####->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
session_name2<-paste(analysis,"_",date,"_","Seurat_pt2.rds",sep="")
session_name3<-paste(analysis,"_",date,"_","Seurat_pt3.rds",sep="")
# DEFAULT PARAMS
figout=TRUE
dimreduce_method="umap"
minpct=0.1#min percentage of cells to have expression of  genes in both groups for DE comparison (here clustering)
logfcthreshold=1# min foldchange (log10 scale) of group2/group1 DE
#
# SL DEFAULT PARAMS
#minpct=0.0125
#logfcthreshold=0.1
#

bm_file1<-"/Users/ngoodacre/Desktop/pipeline_protocols/GEX_Seurat/3Tdata/DFCI.1/ref/combined_celltype_markers.txt"
bm_file2<-"/Users/ngoodacre/Desktop/pipeline_protocols/GEX_Seurat/3Tdata/DFCI.1/ref/SL_additional_markers.txt"
figs_folder_clustered<-paste(figs_folder,"clustered",sep="/")
figs_folder_signs<-paste(figs_folder,"celltype",sep="/")
celltype_signature_ref<-paste(ref_folder,"20200714_cell_signatures_DFCI.1_10x.csv",sep="/")
mref_file<-paste(ref_folder,'20200714_combined_markers.txt',sep="/")
mref<-read.delim(mref_file,sep='\n')[[1]]


r<-RUN_PART_3(session_name2,session_name3,figout,dimreduce_method,minpct,logfcthreshold,bm_file1,bm_file2)
scrna.de<-r[0]
cluster_markers_list<-r[1]
#####->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
##





#########################################################################################################################
#########################################################################################################################
###################################      PART IV -> REPERTOIRE ANALYSIS     #############################################
#########################################################################################################################
#########################################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


RUN_PART_4<-function(session_name3,session_name4,repdat_file){
  ####################################################################################################
  #### LOAD THE DATA FROM PART 3
  scrna.rep=readRDS(paste(sessions_folder,session_name3,sep="/"))#Usually would call this scrna.de, but oh well
  repdat<-read.csv(repdat_file,header=TRUE)
  setup_figs_outs_repertoire(scrna.rep,figs_folder)
  figs_folder_repertoire<-paste(figs_folder,'repertoire',sep='/')
  plot_all_repertoire_info(scrna.rep,repdat,figs_folder_repertoire)
  
  
  #plot_clone_population(vdj,figs_folder_repertoire)

}


##
#####->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
session_name3<-paste(analysis,"_",date,"_","Seurat_pt3.RDS",sep="")
session_name4<-paste(analysis,"_",date,"_","Seurat_pt4.RDS",sep="")
repdat_file="/Users/ngoodacre/Desktop/pipeline_protocols/GEX_Seurat/3Tdata/DFCI.1/ref/20200715_DFCI_repertoire_metadata_Seurat.csv"
RUN_PART_4(session_name3,session_name4,repdat_file)
#####->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
##

