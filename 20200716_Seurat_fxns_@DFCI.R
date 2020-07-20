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
setup_project_env<-function(project_folder,h5_datafolder){
  print(paste("Initializing project ",project_folder,sep=""))
  out_folder=paste(project_folder,"Out",sep="/")
  sessions_folder=paste(project_folder,"Sessions",sep="/")
  figs_folder=paste(out_folder,"Figs",sep="/")
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
  # Detecting organization of data (patient, cell sort, timepoint, multiple files / sample)
  h5_files<-list.files(h5_datafolder,pattern='*.h5')
  num_files<-length(h5_files)
  sample_h5<-h5_files[1]
  s_sample_h5=strsplit(strsplit(sample_h5,"_filtered_feature_bc_matrix.h5")[1][[1]],"_")
  ngroup_levels<-length(s_sample_h5[[1]])
  r<-list()
  r[[1]]<-num_files
  r[[2]]<-ngroup_levels
  return(r)
}





####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
#######################################       PART I        ########################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


####################################################################################################
####################################################################################################
#### SET UP THE SEURAT OBJECT
# Load the raw (non-normalized data).
#scrna.data <- Read10X_h5(data.dir = getwd())
setup_seurat_objects<-function(cwd,h5files,badfiles,ngroups,analysis){
  f=1
  samples<-c()
  seurat_objects=c()
  for (h5file in h5files){
    if(h5file%in%badfiles){
      next
    }
    print(cat("Reading ",h5file))
    s_h5file=strsplit(h5file,"_")[[1]]
    patient=s_h5file[1]
    sort=s_h5file[2]
    if(ngroups==3){
      timepoint=s_h5file[3]
      sample=paste(patient,sort,timepoint,sep="_")
    }
    else{
      sample=paste(patient,sort,sep="_")
    }
    samples=c(samples,sample)
    scrna.data <- Read10X_h5(paste(getwd(),h5file,sep="/"))
    scrna <- CreateSeuratObject(counts = scrna.data, project = sample, min.cells = 3, min.features = 200)
    ncells=ncol(scrna)
    scrna$sample=rep(c(sample),ncells)
    scrna$patient=rep(c(patient),ncells)
    scrna$sort=rep(c(sort),ncells)
    seurat_objects=c(seurat_objects,scrna)
    f=f+1
    seurat_objects[f]=seurat_objects[f][[1]]
  }
  r<-list()
  r[[1]]<-seurat_objects
  r[[2]]<-samples
  return(r)
}


####################################################################################################
#### STANDARD PRE PROCESSING AND NORMALIZATION WORKFLOW
#### This comes before the data ingest and object set-up, because the data have to be QC'ed by sample
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
QC_norm <- function(scrna,figs_folder,name,date){
  print("Performing QC and normalization")
  scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^MT-")
  
  # Show QC metrics for the first 5 cells
  #head(scrna@meta.data, 5)
  
  # Visualize QC metrics as a violin plot
  #VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  
  figs_folder_qcnorm=paste(figs_folder,'QC_norm',sep="/")
  
  if(dir.exists(figs_folder_qcnorm)!=TRUE){
    print(paste("Creating folder ",figs_folder_qcnorm,sep="/"))
    dir.create(figs_folder_qcnorm)
  }
  
  
  pdf(paste(figs_folder_qcnorm,"/",date,"_",name,"_","Preprocessing.pdf",sep=""))
  plot1<-VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(plot1)
  plot2<-FeatureScatter(scrna, feature1 = "nCount_RNA", feature2 = "percent.mt")
  print(plot2)
  plot3<-FeatureScatter(scrna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot3)
  dev.off()
  
  ### implement filtering based on minimum total cell count and percent mitochondrial DNA (dead cell removal)
  #scrna <- subset(scrna, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  ### implement filtering to remove empty cells, multiplet cells (heuristic from Monocle Vignette at: http://cole-trapnell-lab.github.io/monocle-release/docs/)
  scrna.ncounts=scrna$nCount_RNA
  upper_bound <- 10^(mean(log10(scrna.ncounts)) + 2*sd(log10(scrna.ncounts)))# define multiplets
  lower_bound <- 10^(mean(log10(scrna.ncounts)) - 2*sd(log10(scrna.ncounts)))# define empty cells
  scrna=scrna[,scrna.ncounts>lower_bound&scrna.ncounts<upper_bound]
  
  #### NORMALIZING THE DATA
  scrna <- NormalizeData(scrna, normalization.method = "LogNormalize", scale.factor = 10000) #specify defaults
  #scrna <- NormalizeData(scrna) #same
  #scrna <- SCTransform(scrna, vars.to.regress = "percent.mt", verbose = FALSE)# store mitochondrial percentage in object meta data
  #pbmc <- SCTransform(pbmc)# store mitochondrial percentage in object meta data
  scrna <- PercentageFeatureSet(scrna, pattern = "^MT-", col.name = "percent.mt")
  
  ####################################################################################################
  
  return(scrna)
}


####################################################################################################
####################################################################################################
#### FIND WHICH REF/QUERY GENES ARE MISSING IN ANY SEURAT OBJECTS
find_missing_genes<-function(seurat_object_list,refgenes){
  missing<-c()
  for(refgene in refgenes){
    for(seurat_object in seurat_object_list){
      genes<-rownames(seurat_object)
      if(refgene%in%genes==FALSE){
        missing<-c(missing,refgene)
      }
    }
  }
  missing<-unique(missing)
  return(missing)
}


####################################################################################################
####################################################################################################
#### UPDATE ALL SEURAT OBJECTS IN LIST TO CONTAIN REF/PROVIDED LIST (WITH ZEROS)
patch_gene_rows<-function(seurat_object_list,refgenes,analysis){
  new_seurat_object_list=list()
  f=1
  for(seurat_object in seurat_object_list){
    gex<-GetAssayData(seurat_object)
    gex_copy<-gex
    cell_ids<-colnames(gex)
    ncol<-dim(gex)[2]
    hasgenes<-rownames(gex)
    missing_genes<-setdiff(refgenes,hasgenes)
    for(mgene in missing_genes){
      row<-rep(0,ncol)
      names(row)<-cell_ids
      rownames1<-rownames(gex_copy)
      gex_copy<-rbind(gex_copy,row)
      rownames2<-c(rownames1,mgene)
      rownames(gex_copy)<-rownames2
    }
    new_seurat_object<-CreateSeuratObject(counts=gex_copy,project=seurat_object$sample)
    new_seurat_object$sample<-seurat_object$sample
    new_seurat_object$patient<-seurat_object$patient
    new_seurat_object$sort<-seurat_object$sort
    new_seurat_object_list=c(new_seurat_object_list,new_seurat_object)
    new_seurat_object_list[f]<-new_seurat_object_list[f][[1]]
  }
  return(new_seurat_object_list)
}



####################################################################################################
####################################################################################################
#### MERGE THE LOADED DATA
merge_seurat_objects<-function(seurat_objects,samples,ngroup_levels){
  scrna.merged=merge(seurat_objects[[1]], y =c(seurat_objects[[2]],seurat_objects[[3]],seurat_objects[[4]],seurat_objects[[5]],
                                              seurat_objects[[6]],seurat_objects[[7]],seurat_objects[[8]],seurat_objects[[9]],seurat_objects[[10]],
                                              seurat_objects[[11]],seurat_objects[[12]],seurat_objects[[13]],seurat_objects[[14]],seurat_objects[[15]],
                                              seurat_objects[[16]],seurat_objects[[17]],seurat_objects[[18]],seurat_objects[[19]],seurat_objects[[20]],
                                              seurat_objects[[21]],seurat_objects[[22]],seurat_objects[[23]],seurat_objects[[24]],seurat_objects[[25]],
                                              seurat_objects[[26]],seurat_objects[[27]],seurat_objects[[28]],seurat_objects[[29]],seurat_objects[[30]],
                                              seurat_objects[[31]]),add.cell.ids = samples, project = "scrna")
  #### Merge other groups / indiv dataset attributes
  scrna.merged<-add_groups(scrna.merged,ngroup_levels,FALSE,FALSE)
  return(scrna.merged)
}

merge_seurat_objects_DFCI.1<-function(seurat_objects,samples,analysis,ngroup_levels){
  scrna.merged=merge(seurat_objects[[1]], y =c(seurat_objects[[2]],seurat_objects[[3]],seurat_objects[[4]],seurat_objects[[5]],
                                               seurat_objects[[6]],seurat_objects[[7]],seurat_objects[[8]],seurat_objects[[9]],seurat_objects[[10]],
                                               seurat_objects[[11]],seurat_objects[[12]],seurat_objects[[13]],seurat_objects[[14]],seurat_objects[[15]],
                                               seurat_objects[[16]],seurat_objects[[17]],seurat_objects[[18]],seurat_objects[[19]],seurat_objects[[20]],
                                               seurat_objects[[21]],seurat_objects[[22]],seurat_objects[[23]],seurat_objects[[24]],seurat_objects[[25]],
                                               seurat_objects[[26]],seurat_objects[[27]],seurat_objects[[28]],seurat_objects[[29]],seurat_objects[[30]],
                                               seurat_objects[[31]]),add.cell.ids = samples, project = analysis)
  #### Merge other groups / indiv dataset attributes
  scrna.merged<-add_groups(scrna.merged,ngroup_levels,FALSE,FALSE)
  return(scrna.merged)
}


####################################################################################################
####################################################################################################
#### RE-SPLIT MERGED DATA - by PATIENT
split_merged_seurat_object_bypatient<-function(scrna.merged,splitby,figs_folder,date){
  scrna.bypatient<-SplitObject(scrna.merged,split.by=splitby)
  f=1
  for(scrna.patient in scrna.bypatient){
    patient=names(scrna.bypatient)[f]
    scrna.bypatient[patient]<-QC_norm(scrna.patient,figs_folder,patient,date)
    f=f+1
  }
  return(scrna.bypatient)
}


####################################################################################################
####################################################################################################
#### INTEGRATION AND LABEL TRANSFER
## SELECT FEATURES FOR DOWNSTREAM INTEGRATION
integrate_seurat_bypatient<-function(scrna.bypatient,nfeatures,refgenes,npcas,ngroup_levels){
  ##### DEFAULT PARAMS BELOW
  #nfeatures=2000
  #npcas=30
  scrna.features <- SelectIntegrationFeatures(object.list = scrna.bypatient, nfeatures=nfeatures)
  missing_refgenes<-setdiff(refgenes,scrna.features)
  scrna.features<-c(scrna.features,missing_refgenes)
  #### must find k.filter as min cell # across datasets, else default neighbor filter size of 200 will crash anchoors
  k.filter<-min(c(min(sapply(scrna.bypatient,ncol)),200))
  scrna.anchors <- FindIntegrationAnchors(object.list = scrna.bypatient, anchor.features=scrna.features, k.filter=k.filter, dims=1:30)
  #scrna.anchors <- FindIntegrationAnchors(object.list = seurat_objects_list, anchor.features=vargenes, dims=1:30)
  scrna.integrated <- IntegrateData(anchorset = scrna.anchors, dims=1:npcas)
  #### Merge other groups / indiv dataset attributes
  scrna.integrated<-add_groups(scrna.integrated,ngroup_levels,FALSE,FALSE)
  #### Find highly-variable genes again
  scrna.integrated <- FindVariableFeatures(scrna.integrated, selection.method = "vst", nfeatures = nfeatures)
  return(scrna.integrated)
  
  #scrna.list <- PrepSCTIntegration(object.list = seurat_objects_list, anchor.features = scrna.features, ## This is for SCTransform normalization
  # verbose = FALSE)
  ## IDENTIFY ANCHORS AND INTEGRATE DATASETS
  #scrna.anchors <- FindIntegrationAnchors(object.list = seurat_objects_list, normalization.method = "SCT", ## This is for SCTransform normalization
  #                                         anchor.features = scrna.features, verbose = FALSE)
  #scrna.integrated <- IntegrateData(anchorset = scrna.anchors, normalization.method = "SCT", ## This is for SCTransform normalization
  #verbose = FALSE)
}



  

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
#######################################       PART II       ########################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


####################################################################################################
#### SET UP OUTPUT FIGURES FOLDER STRUCTURE - BIOMARKERS
#Can augment this to include other cell group types like timepoint
setup_figs_outs<-function(scrna.clustered,figs_folder){
  figs_folder_global=paste(figs_folder,"global",sep="/")
  if(dir.exists(figs_folder_global)!=TRUE){
    dir.create(figs_folder_global)
    print(paste("Creating ",figs_folder_global,sep=""))
  }
  sorts<-unique(scrna.clustered$timepoint)
  for(sort in sorts){
    figs_folder_sort=paste(figs_folder,sort,sep="/")
    if(dir.exists(figs_folder_sort)!=TRUE){
      dir.create(figs_folder_sort)
      print(paste("Creating ",figs_folder_sort,sep=""))
    }
  }
  patients<-unique(scrna.clustered$patient)
  for(patient in patients){
    figs_folder_patient<-paste(figs_folder,patient,sep="/")
    if(dir.exists(figs_folder_patient)!=TRUE){
      dir.create(figs_folder_patient)
      print(paste("Creating ",figs_folder_patient,sep=""))
    }
    for(sort  in sorts){
      figs_folder_patient_sort<-paste(figs_folder_patient,sort,sep="/")
      if(dir.exists(figs_folder_patient_sort)!=TRUE){
        dir.create(figs_folder_patient_sort)
        print(paste("Creating ",figs_folder_patient_sort,sep=""))
      }
    }
  }
}


####################################################################################################
#### SET UP OUTPUT FIGURES FOLDER STRUCTURE - CELL TYPE / COMBO BIOMARKERS
setup_figs_outs_celltypes<-function(signs,figs_folder){
  figs_folder_CD4<-paste(figs_folder,"CD4",sep="/")
  if(dir.exists(figs_folder_CD4)!=TRUE){
    dir.create(figs_folder_CD4)
    print(paste("Creating ",figs_folder_CD4,sep=""))
  }
  figs_folder_CD8<-paste(figs_folder,"CD8",sep="/")
  if(dir.exists(figs_folder_CD8)!=TRUE){
    dir.create(figs_folder_CD8)
    print(paste("Creating ",figs_folder_CD8,sep=""))
  }
}


####################################################################################################
#### SET UP OUTPUT FIGURES FOLDER STRUCTURE - CLUSTERED OBJECT
setup_figs_outs_clustered<-function(scrna.clustered,figs_folder){
  figs_folder_clustered=paste(figs_folder,"clustered",sep="/")
  if(dir.exists(figs_folder_clustered)!=TRUE){
    dir.create(figs_folder_clustered)
    print(paste("Creating ",figs_folder_clustered,sep=""))
  }
  
  cluster_ids=levels(scrna.clustered$seurat_clusters)
  for(cluster_id in cluster_ids){
    figs_folder_clustered_clusterid=paste(figs_folder_clustered,cluster_id,sep="/")
    if(dir.exists(figs_folder_clustered_clusterid)!=TRUE){
      dir.create(figs_folder_clustered_clusterid)
      print(paste("Creating ",figs_folder_clustered_clusterid,sep=""))
    }
  }
}


explore_dimensionality<-function(scrna.integrated,npcas,figs_folder_global,date){
  ####################################################################################################
  #### DETERMINING THE "DIMENSIONALITY" OF THE DATASET
  # NOTE: This process can take a long time for big datasets, comment out for expediency. More
  # approximate techniques such as those implemented in ElbowPlot() can be used to reduce
  # computation time
  scrna.integrated <- RunPCA(scrna.integrated, features = VariableFeatures(object = scrna.integrated), npcs=npcas)####! SL uses 30, but could be determined by the JackStraw Plot.. only then would have to initially run on a much larger number
  scrna.integrated<- RunUMAP(scrna.integrated, dims = 1:30)
  scrna.integrated <- JackStraw(scrna.integrated, num.replicate = 100)
  scrna.integrated <- ScoreJackStraw(scrna.integrated, dims = 1:20)
  pdf(paste(figs_folder_global,paste(date,'_JackStraw+Elbow.pdf',sep=""),sep="/"))
  print(JackStrawPlot(scrna.integrated, dims = 1:20))
  print(ElbowPlot(scrna.integrated))
  dev.off()
  ####################################################################################################
  #### PERFORMING LINEAR DIMENSIONAL REDUCTION
  # Examine and visualize PCA results a few different ways
  print(scrna.integrated[["pca"]], dims = 1:5, nfeatures = 5)
  pdf(paste(figs_folder_global,paste(date,'_PCA_postreduction.pdf',sep=""),sep="/"))
  VizDimLoadings(scrna.integrated, dims = 1:2, reduction = "pca")
  DimPlot(scrna.integrated, reduction = "pca")
  DimHeatmap(scrna.integrated, dims = 1, cells = 500, balanced = TRUE)
  DimHeatmap(scrna.integrated, dims = 2, cells = 500, balanced = TRUE)
  DimHeatmap(scrna.integrated, dims = 3, cells = 500, balanced = TRUE)
  DimHeatmap(scrna.integrated, dims = 1:15, cells = 500, balanced = TRUE)
  dev.off()
  ####################################################################################################
  return(scrna.integrated)
  ####! Some  automated way of outputting # PCs, i.e. explaining 95% of variability or encompassing all variable genes
  ####################################################################################################
  
}


####################################################################################################
#### CLUSTER THE CELLS
cluster_cells<-function(scrna.integrated,npcas,res){
  # DEFAULT PARAMS
  # npcas=30
  # res=0.5
  scrna.integrated <- FindNeighbors(scrna.integrated, dims = 1:npcas)
  scrna.integrated <- FindClusters(scrna.integrated, resolution = res)
  return(scrna.integrated)
}


####################################################################################################
#### GENERATE GLOBAL UMAP PLOTS 
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn') # or is it 'umap'
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
plot_global_umaps<-function(scrna.clustered,figs_folder_global,groups,date){

  for(group in  groups){
    pdf(paste(figs_folder_global,paste(date,'_UMAP_by=',group,'.pdf',sep=""),sep="/"),width=11,height=8.5)
    ncolors<-length(unlist(unique(scrna.clustered[[group]])))
    if(ncolors<=2){
      colors<-c('#FF0000','#0000FF')[1:ncolors]
    }
    else{
      colors<-DiscretePalette(ncolors,'polychrome')
    }
    embeds = Embeddings(scrna.clustered[[dimreduce_method]])
    minx<-min(embeds[,1])
    maxx=max(embeds[,1])
    miny<-min(embeds[,2])
    maxy<-max(embeds[,2])
    p<-DimPlot(scrna.clustered,group.by=group,label=TRUE,repel=TRUE,pt.size=0.25,label.size=3,cols=colors)+xlim(xmin,xmax)+ylim(ymin,ymax)
    p<-p+coord_fixed(ratio = 1)+ theme(legend.text=element_text(size=12))
    print(p)
    dev.off()
  }
    
    #plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
    #                                                                  override.aes = list(size = 3)))

}
####################################################################################################





####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
######################################       PART III       ########################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


#########################################################################################################################
#### VISUALIZING IMMUNE BIOMARKER EXPRESSION LEVELS
visualize_biomarkers<-function(scrna.clustered,figs_folder,dimreduce_method,bm_file1,bm_file2,date){
  # DEFAULT PARAM 
  # bm_file="/Users/ngoodacre/Desktop/pipeline_protocols/GEX_Seurat/3Tdata/MERCK.1/SL_analysis_2019/marker_genes.txt" # Shane's / 3T's  list of markers  #some numbering is in col1
  biomarkers1=read.table(bm_file1,sep="\t")[,1]
  biomarkers2=read.table(bm_file2,sep="\t")[,1]
  combined_biomarkers<-c(biomarkers1,biomarkers2)
  for(biomarker in combined_biomarkers){
    print(paste("Plotting","global",biomarker,"expression",sep=" "))
    figs_folder_global<-paste(figs_folder,"global",sep="/")
    plot_nonzeromarker_embeds(scrna.clustered,dimreduce_method,biomarker,TRUE,1.0,figs_folder_global,paste("global","im","alpha=1.0",sep="_"),date)
    plot_nonzeromarker_embeds(scrna.clustered,dimreduce_method,biomarker,TRUE,0.25,figs_folder_global,paste("global","im","alpha=0.25",sep="_"),date)
    #sorts<-unique(scrna.clustered$sort)
    sorts<-unique(scrna.clustered$timepoint)
    for(sort in sorts){
      print(paste("Plotting",sort,biomarker,"expression",sep=" "))
      figs_folder_sort<-paste(figs_folder,sort,sep="/")
      scrna.clustered.sort<-get_groups(scrna.clustered,'timepoint',sort)
      plot_nonzeromarker_embeds(scrna.clustered.sort,dimreduce_method,biomarker,TRUE,1.0,figs_folder_sort,paste(sort,"im","alpha=1.0",sep="_"),date)
      plot_nonzeromarker_embeds(scrna.clustered.sort,dimreduce_method,biomarker,TRUE,0.25,figs_folder_sort,paste(sort,"im","alpha=0.25",sep="_"),date)
    }
    patients<-unique(scrna.clustered$patient)
    for(patient in patients){
      scrna.clustered.patient<-get_groups(scrna.clustered,'patient',patient)
      for(sort in sorts){
        if(sort%in%unique(scrna.clustered.patient$timepoint)){
          print(paste("Plotting",patient,sort,biomarker,"expression",sep=" "))
          figs_folder_patient_sort<-paste(figs_folder,patient,sort,sep="/")
          scrna.clustered.patient.sort<-get_groups(scrna.clustered.patient,'timepoint',sort)
          plot_nonzeromarker_embeds(scrna.clustered.patient.sort,dimreduce_method,biomarker,TRUE,1.0,figs_folder_patient_sort,paste(patient,sort,"im","alpha=1.0",sep="_"),date)
          plot_nonzeromarker_embeds(scrna.clustered.patient.sort,dimreduce_method,biomarker,TRUE,0.25,figs_folder_patient_sort,paste(patient,sort,"im","alpha=0.25",sep="_"),date)
          }
        else{
          print(paste("Patient",patient,"does not have cell sort",sort,sep=" "))
        }
      }
    }
  }
  for(groupby in c('orig.ident','patient','sort','timepoint','patient_timepoint')){#add timepoint,patient_timepoint
    fig_filename<-fig_filename_maker(figs_folder_global,date,paste('groupby=',groupby,sep=""),'im1_Dotplot')
    dotplot_groups(scrna.clustered,groupby,biomarkers1,fig_filename)
    fig_filename<-fig_filename_maker(figs_folder_global,date,paste('groupby=',groupby,sep=""),'im2_Dotplot')
    dotplot_groups(scrna.clustered,groupby,biomarkers2,fig_filename)
  }
}

#########################################################################################################################
#### VISUALIZING IMMUNE BIOMARKER EXPRESSION LEVELS
visualize_biomarkers_dotplot<-function(scrna.clustered,figs_folder,dimreduce_method,bm_file1,bm_file2,date){
  # DEFAULT PARAM 
  # bm_file="/Users/ngoodacre/Desktop/pipeline_protocols/GEX_Seurat/3Tdata/MERCK.1/SL_analysis_2019/marker_genes.txt" # Shane's / 3T's  list of markers  #some numbering is in col1
  biomarkers1=read.table(bm_file1,sep="\t")[,1]
  biomarkers2=read.table(bm_file2,sep="\t")[,1]
  combined_biomarkers<-c(biomarkers1,biomarkers2)
  for(groupby in c('orig.ident','patient','sort','timepoint','patient_timepoint')){#add timepoint,patient_timepoint
    fig_filename<-fig_filename_maker(figs_folder_global,date,paste('groupby=',groupby,sep=""),'im1_Dotplot')
    dotplot_groups(scrna.clustered,groupby,biomarkers1,fig_filename)
    fig_filename<-fig_filename_maker(figs_folder_global,date,paste('groupby=',groupby,sep=""),'im2_Dotplot')
    dotplot_groups(scrna.clustered,groupby,biomarkers2,fig_filename)
  }
}


#########################################################################################################################
#### FINDING DIFFERENTIALLY EXPRESSED CLUSTER BIOMARKERS
#### Find  all cluster markers
find_clustermarkers<-function(scrna.clustered,figout,dimreduce_method,minpct,logfcthreshold,figs_folder_clustered,date){
  cluster_ids=levels(scrna.clustered$seurat_clusters)
  cluster_markers_list=list()
  for(cluster_id in cluster_ids){
    cluster.markers <- FindMarkers(scrna.clustered, ident.1 = cluster_id, min.pct = minpct , logfc.threshold = logfcthreshold)
    cluster_markers_list[cluster_id]<-cluster.markers
    figs_folder_clustered_clusterid=paste(figs_folder_clustered,cluster_id,sep="/")
    cluster.markers.genenames<-rownames(cluster.markers)
    if(figout){
      for(clustermarker in  cluster.markers.genenames){
        plot_nonzeromarker_embeds(scrna.clustered,dimreduce_method,clustermarker,TRUE,1.0,figs_folder_clustered_clusterid,paste(cluster_id,"cm","alpha=1.0",sep="_"),date)
        plot_nonzeromarker_embeds(scrna.clustered,dimreduce_method,clustermarker,TRUE,0.25,figs_folder_clustered_clusterid,paste(cluster_id,"cm","alpha=0.25",sep="_"),date)
      }
    }
    for(groupby in c('orig.ident','patient','sort')){#add timepoint,patient_timepoint
      fig_filename<-fig_filename_maker(figs_folder_clustered_clusterid,date,paste('groupby=',groupby,sep=""),paste("cluster=",cluster_id,'_Dotplot',sep=""))
      dotplot_groups(scrna.clustered,groupby,cluster.markers.genenames,fig_filename)
    }
  }
  return(c(scrna.clustered,cluster_markers_list))
}


#########################################################################################################################
#### DEFINING ANDS VISUALIZING IMMUNE CELL POPULATIONS (CD8 naive, effector, memory, (pre)dysfunctional ; 
#### CD4 naive, helper, memory, resident -


parse_celltype_signatures<-function(celltype_signature_ref){
  sign_table<-read.csv(celltype_signature_ref,header=TRUE)
  cd_4_8<-sign_table[,'CD4_or_CD8']
  refs<-sign_table[,'dataset']
  sign_table<-sign_table[cd_4_8=='CD4' & refs=='Geginat et al'|cd_4_8=='CD8' & refs=='3T_unified',]
  signs<-list()
  for(r in 1:nrow(sign_table)){
    row<-sign_table[r,]
    cd_4_8<-toString(row['CD4_or_CD8'])
    pop=toString(row['population'])
    #subpop=row['subpopulation']#currently under development (20200714)
    pos_markers=strsplit(toString(row['positive_gene_signature'])[[1]],", ")
    neg_markers=strsplit(toString(row['negative_gene_signature'])[[1]],", ")
    signs[[cd_4_8]][[pop]][['+']]<-pos_markers
    signs[[cd_4_8]][[pop]][['-']]<-neg_markers
  
  }
  return(signs)
}

define_CD4_CD8<-function(scrna.clustered){
  gex<-GetAssayData(scrna.clustered)
  cd8a<-names(get_nonzero_marker_cells(gex,'CD8A',TRUE))#Taking union here
  cd8b<-names(get_nonzero_marker_cells(gex,'CD8B',TRUE))
  cd8<-list()
  cd8[['CD8A']]=cd8a
  cd8[['CD8B']]=cd8b
  cd8<-c(cd8a,cd8b)
  cd4<-names(get_nonzero_marker_cells(gex,'CD4',TRUE))
  allcells=colnames(gex)
  cd4<-allcells%in%cd4
  cd8<-allcells%in%cd8
  scrna.clustered$CD4<-cd4
  scrna.clustered$CD8<-cd8
  return(scrna.clustered)
}


define_cell_types<-function(scrna.clustered,signs){
  scrna.clustered<-define_CD4_CD8(scrna.clustered)
  allcells<-colnames(scrna.clustered)
  cell_types<-list()
  for(cd_4_8 in names(signs)){
    gex<-GetAssayData(scrna.clustered)
    gex_cd<-gex[,which(scrna.clustered[[cd_4_8]][[1]]==TRUE)]
    signs_cd<-signs[[cd_4_8]]
    for(pop in names(signs_cd)){
      signs_pop<-signs_cd[[pop]]
      signs_pop.pos<-signs_pop[['+']][[1]]
      signs_pop.neg<-signs_pop[['-']][[1]]
      pop_cells<-list()
      for(marker in signs_pop.pos){
        marker_cells<-names(get_nonzero_marker_cells(gex_cd,marker,TRUE))#last param indicates whether to positively select
        pop_cells[[marker]]<-marker_cells
      }
      for(marker in signs_pop.neg){
        marker_cells<-names(get_nonzero_marker_cells(gex_cd,marker,FALSE))#last param indicates whether to positively select
        pop_cells[[marker]]<-marker_cells
      }
      pop_cells<-Reduce(intersect,pop_cells)
      cell_type_membership<-allcells%in%pop_cells#dataset-wide vector, boolean for defining cell type/pop membership
      names(cell_type_membership)<-allcells
      cell_types[[cd_4_8]][[pop]]<-cell_type_membership
      pop_fullname<-paste(cd_4_8,pop,sep="_")
      scrna.clustered[[pop_fullname]]=cell_type_membership
    }
  }
  r<-list()
  r[[1]]<-scrna.clustered
  r[[2]]<-cell_types
  return(r)
}


plot_cell_types<-function(scrna.clustered,signs,filter_cell_ids,filter_name,figs_folder_signs,date){
  for(cd_4_8 in names(signs)){
    # Plot CD4 and CD8 cells
    membership<-scrna.clustered[[cd_4_8]][[1]]
    member_cell_ids<-colnames(scrna.clustered)[membership==TRUE]
    membership<-ifelse(membership==TRUE,cd_4_8,'other')#one-hot encode membership to cell population
    if(length(filter_cell_ids)>0){
      cell_ids_list<-list()
      cell_ids_list[[1]]=member_cell_ids
      cell_ids_list[[2]]=filter_cell_ids
      member_cell_ids<-Reduce(intersect,cell_ids_list)
    }
    scrna.clustered.copy<-scrna.clustered
    scrna.clustered.copy[[cd_4_8]]=membership
    print(paste("Plotting cell type ",cd_4_8,sep=" "))
    umap_bycells(scrna.clustered.copy,date,member_cell_ids,filter_name,cd_4_8, c('red','grey'),figs_folder_signs)
    
    signs_cd<-signs[[cd_4_8]]
    for(pop in names(signs_cd)){
      
      pop_fullname<-paste(cd_4_8,pop,sep="_")
      print(paste("Plotting cell type ",pop_fullname,sep=" "))
      membership<-scrna.clustered[[pop_fullname]][[1]]
      member_cell_ids<-colnames(scrna.clustered)[membership==TRUE]
      if(length(filter_cell_ids)>0){
        cell_ids_list<-list()
        cell_ids_list[[1]]=member_cell_ids
        cell_ids_list[[2]]=filter_cell_ids
        member_cell_ids<-Reduce(intersect,cell_ids_list)
      }
      if(length(member_cell_ids)<=1){
        print(paste("Insufficient size for",pop_fullname,"population - skipping",sep=" "))
        next
      }


      membership<-ifelse(membership==TRUE,pop_fullname,'other')#one-hot encode membership to cell population
      scrna.clustered.copy<-scrna.clustered
      scrna.clustered.copy[[pop_fullname]]=membership
      umap_bycells(scrna.clustered.copy,date,member_cell_ids,filter_name,pop_fullname,c(),figs_folder_signs)
    }
  }
}


plot_cell_types__main<-function(scrna.clustered,signs,figs_folder,date){
  print(paste("Plotting","global","cell types",sep=" "))
  figs_folder_global<-paste(figs_folder,"global",sep="/")
  plot_cell_types(scrna.clustered,signs,c(),"global",figs_folder_global,date)
  timepoints<-unique(scrna.clustered$timepoint)
  
  for(timepoint in timepoints){
    print(paste("Plotting",timepoint,"cell types",sep=" "))
    figs_folder_timepoint<-paste(figs_folder,timepoint,sep="/")
    scrna.clustered.timepoint<-scrna.clustered[,scrna.clustered$timepoint==timepoint]
    filter_cells<-colnames(scrna.clustered.timepoint)
    plot_cell_types(scrna.clustered,signs,filter_cells,timepoint,figs_folder_timepoint,date)
  }
  
  patients<-unique(scrna.clustered$patient)
  for(patient in patients){
    scrna.clustered.patient<-get_groups(scrna.clustered,'patient',patient)
    for(timepoint in timepoints){
        if(timepoint%in%unique(scrna.clustered.patient$timepoint)){
          patient_timepoint<-paste(patient,timepoint,sep="_")
          print(paste("Plotting",patient,timepoint,"cell types",sep=" "))
          figs_folder_patient_timepoint<-paste(figs_folder,patient,timepoint,sep="/")
          scrna.clustered.patient.timepoint<-get_groups(scrna.clustered.patient,'timepoint',timepoint)
          filter_cells<-colnames(scrna.clustered.patient.timepoint)
          plot_cell_types(scrna.clustered,signs,filter_cells,patient_timepoint,figs_folder_patient_timepoint,date)
        }
        else{
          print(paste("Patient",patient,"does not have timepoint",timepoint,sep=" "))
        }
      }
    }#
}

  
  



####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
##########################       PART IV REPERTOIRE ANALYSIS       #################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################





####################################################################################################
#### SET UP OUTPUT FIGURES FOLDER STRUCTURE - CLUSTERED OBJECT
setup_figs_outs_repertoire<-function(scrna.rep,figs_folder){
  figs_folder_repertoire=paste(figs_folder,"repertoire",sep="/")
  if(dir.exists(figs_folder_repertoire)!=TRUE){
    dir.create(figs_folder_repertoire)
    print(paste("Creating ",figs_folder_repertoire,sep=""))
  }
  rep_analyses<-c("expanded","persistent","TCR_clusters")
  for(rep_analysis in rep_analyses){
    figs_folder_repertoire_analysis<-paste(figs_folder_repertoire,rep_analysis,sep="/")
    if(dir.exists(figs_folder_repertoire_analysis)!=TRUE){
      dir.create(figs_folder_repertoire_analysis)
      print(paste("Creating ",figs_folder_repertoire_analysis,sep=""))
    }
  }

}


# Filter repertoire dataframe by patient / timepoint (e.g. persistence) / clustered criteria
filter_repdat_1<-function(repdat,col,filters){
  if(col=='patient'){
    subrepdat<-subset(repdat,patient%in%filters)
  }
  if(col=='presence'){#'presence' = whether persistent, pre, or post
    subrepdat<-subset(repdat,presence%in%filters)
  }
  if(col=='cluster'){
    subrepdat<-subset(repdat,cluster%in%filters)
    }
  return(subrepdat)
}

# Filter repertoire dataframe by clonal expansion criteria
filter_repdat_2<-function(repdat,metric,threshold){
  patients<-repdat[,'patient']
  #pre_clonect<-repdat[,'predose_cellcount']
  #post_clonect<-repdat[,'postdose_cellcount']
  #pre_clonefrac<-repdat[,'predose_cellfrac']
  #post_clonefrac<-repdat[,'postdose_cellfrac']
  #pre_clonerank<-repdat[,'predose_cellrank']
  #post_clonerank<-repdat[,'postdose_cellrank']
  
  if(metric=='predose_cellcount'){
    subrepdat<-subset(repdat,predose_cellcount>threshold)
  }
  if(metric=='postdose_cellcount'){
    subrepdat<-subset(repdat,postdose_cellcount>threshold)
  }
  
  if(metric=='predose_cellfrac'){
    subrepdat<-subset(repdat,predose_cellfrac>threshold)
  }
  if(metric=='postdose_cellfrac'){
    subrepdat<-subset(repdat,postdose_cellfrac>threshold)
  }
  
  if(metric=='predose_rank'){
    subrepdat<-subset(repdat,predose_rank<=threshold)
  }
  if(metric=='postdose_rank'){
    subrepdat<-subset(repdat,postdose_rank<=threshold)
  }
  
  return(subrepdat)

}


# Filter repertoire dataframe by specific clone(s) criteria
filter_repdat_3<-function(repdat,col,filters){
  subrepdat<-subset(repdat,clonotype%in%filters)
  return(subrepdat)
}


merge_repdats<-function(repdat1,repdat2,repdat_ids){
  repdat.merged=merge(scrna.rep[[1]], y =c(seurat_objects[[2]]),add.cell.ids = repdat_ids, project = "rep")
  return(repdat.merged)
}


# Extract VDJ-repertoire cellids (concatenated patient, sort (assumes either CD3 or ALL), timepoint, cellular barcode) for repdat
extract_repdat_cellids<-function(repdat){
  repdat_cellids<-c()
  cellids_col<-repdat[,'seurat_cellids']
  for(cellids in cellids_col){
    s_cellids<-strsplit(cellids," ; ")[[1]]
    repdat_cellids<-c(repdat_cellids,s_cellids)
  }
  return(repdat_cellids)
}


# For looking at intersection between Seurat and VDJ-repertoire barcodes (given potential groups in former)
intersect_repdat_seurat_cellids_OLD<-function(seurat_cellids,rep_cellids){#! obsolete, unecessary, slow
  isect=c()
  for(sbc in seurat_cellids){
    match=FALSE
    if(sbc%in%rep_cellids){
      match=TRUE
    }
    s_sbc=strsplit(sbc,"-")
    if(length(s_sbc)>2){
      short_sbc=paste(s_sbc[1],s_sbc[2],sep="-")
      if(short_sbc%in%rep_cellids){
        match=TRUE
      }
    }
    if(match==TRUE){
      isect<-c(isect,sbc)
    }
  }
  return(isect)
}

intersect_repdat_seurat_cellids<-function(seurat_cellids,rep_cellids){
  idlist<-list()
  idlist[[1]]=seurat_cellids
  idlist[[2]]=rep_cellids
  isect<-Reduce(intersect,idlist)
  return(isect)
}


# Essentially combines the previous two functions
# extract_repdat_cellids + intersect_barcoded_cellids
repdat_to_seurat<-function(scrna.rep,subrepdat){
  seurat_cellids<-colnames(scrna.rep)
  repdat_cellids<-extract_repdat_cellids(subrepdat)
  seurat_cellids_x<-intersect_repdat_seurat_cellids(seurat_cellids,repdat_cellids)
  return(seurat_cellids_x)
}


# Plotting function to make a cell-specific umap for a single repertoire cross-section (binary for criterion)
plot_repertoire_info_binary<-function(scrna.rep,subrepdat,inclass,outclass,classname,date,figs_folder_rep){
  cell_ids<-repdat_to_seurat(scrna.rep,subrepdat)
  membership<-cellids_membership(scrna.rep,cell_ids,inclass,outclass)
  scrna.rep[[classname]]<-membership
  umap_bycells(scrna.rep,date,cell_ids,inclass,classname,c(),figs_folder_rep)
  return(scrna.rep)
}


# Plotting function to make a cell-specific umap for a single repertoire cross-section (multi-class )
plot_repertoire_info_multi<-function(scrna.rep,subrepdat_vect,inclass_vect,outclass,classname,date,figs_folder_rep){
  memberships<-list()
  for(i in 1:length(inclass_vect)){
    inclass<-inclass_vect[i]
    subrepdat<-subrepdat_vect[[i]]
    cell_ids<-repdat_to_seurat(scrna.rep,subrepdat)
    membership<-cellids_membership(scrna.rep,cell_ids,inclass,outclass)
    memberships[[i]]<-membership
  }
  merged_membership<-merge_cellids_membership(memberships,outclass)
  i1<-which(merged_membership=='pre+persistent')
  i2<-which(merged_membership=='post+persistent')
  merged_membership<-replace(merged_membership,i1,rep("persistent",length(i1)))
  merged_membership<-replace(merged_membership,i2,rep("persistent",length(i2)))
  scrna.rep[[classname]]<-merged_membership
  cell_ids<-colnames(scrna.rep)[which(merged_membership!=outclass)]
  colors<-DiscretePalette(4,'alphabet')
  colors<-DiscretePalette(3,'alphabet')
  #colors[4]<-"#eeeeee"
  umap_bycells(scrna.rep,date,cell_ids,classname,classname,colors,figs_folder_rep)
  return(scrna.rep)
}


# Plotting function to make a cell-specific umap for a two repertoire cross-sections (multi-class)
# where one is visualized by color (colorgroup), the other by shape (shapegroup)
#!! also want to add the ability to highlight certain cells and resize
plot_repertoire_info_multi2<-function(scrna.rep,subrepdat_list,inclass_vect,outclass,classname,shapeby,date,figs_folder_rep){
  memberships<-list()
  cell_ids_sub<-c()
  for(i in 1:length(inclass_vect)){
    inclass<-inclass_vect[i]
    subrepdat<-subrepdat_list[[i]]
    cell_ids<-repdat_to_seurat(scrna.rep,subrepdat)
    cell_ids_sub<-c(cell_ids_sub,cell_ids)
    membership<-cellids_membership(scrna.rep,cell_ids,inclass,outclass)
    memberships[[i]]<-membership
  }
  cell_ids_sub<-unique(cell_ids_sub)
  merged_membership<-merge_cellids_membership(memberships,outclass)
  scrna.rep[[classname]]<-merged_membership
  highlight_cellids<-c()
  colors<-DiscretePalette(3,'alphabet')
  umap_bycells2(scrna.rep,date,cell_ids_sub,classname,classname,shapeby,highlight_cellids,colors,figs_folder)
  umap_bycells2(scrna.rep,date,cell_ids_sub,classname,'patient',shapeby,highlight_cellids,c(),figs_folder)
    
  return(scrna.rep)
}


################################################################################################
################################################################################################
# Control function to plot all basic repertoire info
plot_all_repertoire_info(scrna.rep,repdat,figs_folder_rep){
  
  # Plot highly-expanded clones
  #repdat_expanded.post.1<-filter_repdat_2(repdat,'postdose_rank',1)
  #repdat_expanded.post.10<-filter_repdat_2(repdat,'postdose_rank',10)
  #scrna.rep<-plot_repertoire_info_binary(scrna.rep,repdat_expanded.post.rank.1,'expanded_top1','less_expanded','r1',date,figs_folder_rep)
  #scrna.rep<-plot_repertoire_info_binary(scrna.rep,repdat_expanded.post.rank.10,'expanded_top10','less_expanded','r10',date,figs_folder_rep)
  
  # Plot highly-expanded clones
  repdat_expanded.post.10<-filter_repdat_2(repdat,'postdose_cellcount',10)
  repdat_expanded.post.100<-filter_repdat_2(repdat,'postdose_cellcount',100)
  scrna.rep<-plot_repertoire_info_binary(scrna.rep,repdat_expanded.post.10,'expanded_post10','less_expanded','expanded_post10',date,figs_folder_rep)
  scrna.rep<-plot_repertoire_info_binary(scrna.rep,repdat_expanded.post.100,'expanded_post100','less_expanded','expanded_post100',date,figs_folder_rep)
  
  # Plot pre/post/persistent clones
  repdat_timepoint.pre<-filter_repdat_1(repdat,'presence','pre')
  repdat_timepoint.post<-filter_repdat_1(repdat,'presence','post')
  repdat_timepoint.persistent<-filter_repdat_1(repdat,'presence','persistent')
  scrna.rep<-plot_repertoire_info_multi(scrna.rep,c(repdat_timepoint.pre,repdat_timepoint.post,repdat_timepoint.persistent),
                             c('pre','post','persistent'),'unknown','persistence',date,figs_folder_rep)
  
  # Plot persistent clones with their non-persistent cluster members
  repdat_clusters.C1<-filter_repdat_1(repdat,'cluster',c('C1'))
  repdat_clusters.C2<-filter_repdat_1(repdat,'cluster',c('C2'))
  repdat_clusters_list<-list()
  repdat_clusters_list[[1]]<-repdat_clusters.C1
  repdat_clusters_list[[2]]<-repdat_clusters.C2
  scrna.rep<-plot_repertoire_info_multi2(scrna.rep,repdat_clusters_list,c('C1','C2'),'non/minor-cluster','TCRcluster','persistence',date,figs_folder_rep)
  
  # Plot "head nodes" of clusters - the highly-expanded central, persistent clones
  repdat_clusters.C1.headNOVMEL5<-filter_repdat_3(repdat,'clonotype',c('CASSPPQGAQKLVF__CASSRPELDTQYF'))
  repdat_clusters.C2.headNEIMEL1<-filter_repdat_3(repdat,'clonotype',c('CAVRYGQNFVF__CASSLEGDRPQHF'))
  repdat_clusters.C2.headNOVMEL3<-filter_repdat_3(repdat,'clonotype',c('CAVRRGQNFVF__CASSSEGDRPQYF'))
  scrna.rep<-plot_repertoire_info_binary(scrna.rep,repdat_clusters.C1.headNOVMEL5,'C1.headNOVMEL5','other','C1.headNOVMEL5',date,figs_folder_rep)
  scrna.rep<-plot_repertoire_info_binary(scrna.rep,repdat_clusters.C2.headNEIMEL1,'C2.headNEIMEL1','other','C2.headNEIMEL1',date,figs_folder_rep)
  scrna.rep<-plot_repertoire_info_binary(scrna.rep,repdat_clusters.C2.headNOVMEL3,'C2.headNOVMEL3','other','C2.headNOVMEL3',date,figs_folder_rep)
  
  return(scrna.rep)

}



























# Reads in metadata for VDJ repertoire; stores by sample_barcode
get_VDJ_metadata<-function(vdj_file){
  vdj_data=read.csv(vdj_file,header=TRUE)
  vdj=list()
  for(r in 1:nrow(vdj_data)){
    row=vdj_data[r,]
    bc=row['barcode']
    group=row['group']
    rank=as.numeric(row['rank'])
    row['rank']=rank
    cell_id=paste(group,bc,sep="_")
    vdj[cell_id]=list(row)
  }
  return(vdj)
}

# Takes in param / value / fun(param) and on VDJ metadata, returns only barcodes matching those params
filter_cells<-function(vdj,param,val,func){
  vdj2=list()
  cell_ids<-names(vdj)
  r=1
  for(row in vdj){
    val01=row[param]
    print(row['rank'])
    if(func(val01,val)){
      cell_id<-cell_ids[r]
      vdj2[cell_id]=list(row)
    }
    r=r+1
  }
  return(vdj2)
}

# Creates a column vector for a Seurat object with 1=present, 0=absent in filter set(cell_ids)
cell_cols<-function(scrna,cell_ids){
  cell_col<-ifelse(colnames(scrna)%in%cell_ids,1,0)
  return(cell_col)
}


# Takes seurat object as input, returns cells barcodes
get_barcodes<-function(scrna){
  mat=GetAssayData(scrna)
  barcodes=c()
  cell_ids<-rownames(mat)
  for(cell_id in cell_ids){
    s_cell_id=strsplit(cell_id,"_")[[1]]
    barcode=s_cell_id[length(s_cell_id)]
    barcodes=c(barcodes,barcode)
  }
  return(barcodes)
}


ingroup_gex<-function(scrna,cell_ids,marker){
  scrna.2<-scrna[,cell_ids]
  marker_gex<-scrna.1[marker,]
  marker_gex_mean<mean(marker_gex)
  return(marker_gex_mean)
}




plot_clone_population<-function(vdj,figs_folder_repertoire){
  vdj2=filter_cells(vdj,"rank",10,lessthan)
  vdj3=filter_cells(vdj,"rank",100,lessthan)
  vdj4=filter_cells(vdj,"persistent","Y",equals)
  vdj5=filter_cells(vdj,"persistent_expinpost",1,greaterthan)
  vdj6=filter_cells(vdj,"persistent_expinpost",5,greaterthan)
  vdj7=filter_cells(vdj,"tcr-dist_cluster","c1",equals)
  vdj7=filter_cells(vdj,"tcr-dist_cluster","c2",equals)
  ##### ALso need to plot all clusters simultaneously, by color
  pop_names<-c('expandedtop10','expandedtop100','persistent','persistent_expinpost_>1','persistent_expinpost_>5','TCR_cluster1','TCR_cluster2')
  i=1
  for(vdjx in c(vdj2,vdj3,vdj4,vdj5,vdj6,vdj7)){
    pop_name<-pop_names[i]
    plot_clone_population(scrna.de,vdjx,pop_name)
    i<-i+1
    cell_ids<-names(vdjx)
    cell_ids_name<-'expandedtop10'
    for(groupby in c('sort','patient','orig.ident')){
      umap_bycells(scrna,date,cell_ids,cell_ids_name,groupby,c(),figs_folder_repertoire)
    }
  }

}































#########################################################################################################################
#########################################################################################################################
####################################      VISUALIZATION FUNCTIONS    ####################################################
#########################################################################################################################
#########################################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


library(viridis)

# scrna = main seurat object  for analysis (integrated across patients, likely also dim-reduced, clustered)
# dimreduce_method =  e.g. "umap", "tsne"
# marker = gene
plot_nonzeromarker_embeds<-function(scrna,dimreduce_method,marker,filterzero,transparency,figs_folder_X,handle,date){
  gex=GetAssayData(scrna)
  gex_marker=gex[rownames(gex)==marker,]
  if(length(gex_marker)==0){
    return('no_cells_with_marker')
  }
  embeds = Embeddings(scrna[[dimreduce_method]])
  minx<-min(embeds[,1])
  maxx=max(embeds[,1])
  miny<-min(embeds[,2])
  maxy<-max(embeds[,2])
  #embeds = scrna.pt1@reductions$umap
  # Filter out empty cells (if "filterzero"==TRUE)
  if(filterzero){
    lower_bound <- mean(gex_marker) + sd(gex_marker)# define boundary for empty cells (cells with <  1SD ABOVE MEAN)
    embeds=embeds[gex_marker>lower_bound,]
    gex_marker=gex_marker[gex_marker>lower_bound]
  }
  embeds=data.frame(embeds)
  embeds$u1<-embeds[,1]
  embeds$u2<-embeds[,2]

  embeds$gex<-gex_marker#this would be gene-specific
  pdf(paste(figs_folder_X,paste(paste(date,handle,dimreduce_method,marker,sep="_"),'.pdf',sep=""),sep="/"))
  p<-ggplot(embeds,aes(x=u1,y=u2))
  print(p+geom_point(aes(col=gex_marker))+scale_color_viridis(option="C",direction=1,alpha=transparency)+xlim(minx,maxx)+ylim(miny,maxy))
  dev.off()
  
}


# Plot umap with only certain cells highlighted, retaining dimensions of original UMAP
umap_bycells<-function(scrna,date,cell_ids,cell_ids_name,groupby,colors,figs_folder){
  #cols=c('#F8766D','#00BA38')
  
  embeds = Embeddings(scrna[['umap']])
  u1=embeds[,1]
  u2=embeds[,2]
  xmin=min(u1)
  xmax=max(u1)
  ymin=min(u2)
  ymax=max(u2)
  
  scrna_slice<-scrna[,cell_ids]
  pdf(paste(figs_folder,paste(date,'_UMAP_by=',groupby,'_cells=',cell_ids_name,'.pdf',sep=""),sep="/"),width=11,height=8.5)
  if(length(colors)==0){
    ncolors<-length(unique(unlist(scrna_slice[[groupby]])))
    if(ncolors<=2){
      colors<-c('#FF0000','#0000FF')[1:ncolors]
    }
    else{
      colors<-DiscretePalette(ncolors,'polychrome')
    }
  }
  p<-DimPlot(scrna_slice,group.by=groupby,label=TRUE,repel=TRUE,pt.size=1,label.size=3,cols=colors)+xlim(xmin,xmax)+ylim(ymin,ymax)
  p<-p+coord_fixed(ratio = 1)+ theme(legend.text=element_text(size=14))

  print(p)

  #plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
  #                                                                    override.aes = list(size = 3)))
  dev.off()
}


# Plot umap with only certain cells highlighted, retaining dimensions of original UMAP#!! Replace previous with this, post 20200716
umap_bycells2<-function(scrna,date,cell_ids,cell_ids_name,groupby,shapeby,highlight_cellids,colors,figs_folder){
  #cols=c('#F8766D','#00BA38')
  embeds = Embeddings(scrna[['umap']])
  u1=embeds[,1]
  u2=embeds[,2]
  xmin=min(u1)
  xmax=max(u1)
  ymin=min(u2)
  ymax=max(u2)
  
  scrna_slice<-scrna[,cell_ids]
  pdf(paste(figs_folder,paste(date,'_UMAP_by=',groupby,'_cells=',cell_ids_name,'_shapeby=',shapeby,'.pdf',sep=""),sep="/"),width=11,height=8.5)
  if(length(colors)==0){
    ncolors<-length(unique(unlist(scrna_slice[[groupby]])))
    if(ncolors<=2){
      colors<-c('#FF0000','#0000FF')[1:ncolors]
    }
    else{
      colors<-DiscretePalette(ncolors,'polychrome')
    }
  }
  p<-DimPlot(scrna_slice,group.by=groupby,shape.by=shapeby,label=TRUE,repel=TRUE,pt.size=2,label.size=3,cols=colors)+xlim(xmin,xmax)+ylim(ymin,ymax)
  p<-p+coord_fixed(ratio = 1)+ theme(legend.text=element_text(size=14))
  print(p)
  
  #plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
  #                                                                    override.aes = list(size = 3)))
  dev.off()
}


fig_filename_maker<-function(figs_folder,date,tag,analysis){
  fig_filename<-paste(figs_folder,paste(paste(date,tag,analysis,sep="_"),".pdf",sep=""),sep="/")
  return(fig_filename)
}


dotplot_groups<-function(scrna,groupby,markers,fig_filename){
  pdf(fig_filename)
  print(DotPlot(scrna,assay="RNA",group.by=groupby,features=markers,col.min=-10,dot.min=0)+RotatedAxis()+FontSize(y.text=6)+coord_flip())
  dev.off()
}













#########################################################################################################################
#########################################################################################################################
########################################      GENERIC FUNCTIONS    ######################################################
#########################################################################################################################
#########################################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


# Takes seurat object as input, returns cells barcodes
get_barcodes<-function(scrna){
  mat=GetAssayData(scrna)
  barcodes=c()
  for(cell_id in cell_ids){
    s_cell_id=strsplit(cell_id,"_")[[1]]
    barcode=s_cell_id[length(s_cell_id)]
    barcodes=c(barcodes,barcode)
  }
  return(barcodes)
}


# parameter-passable function for "<" operator ; is there another way?
lessthan<-function(x,y){
  return(x<y)
}


# parameter-passable function for ">" operator ; is there another way?
morethan<-function(x,y){
  return(x>y)
}

equals<-function(x,y){
  return(x==y)
}

# Merge other groups / indiv dataset attributes - MERCK DATASET (deprecated)
add_groups_MERCK<-function(scrna.big){
  samples.big=scrna.big$orig.ident
  patients.big=c()
  sorts.big=c()
  for(sample in samples.big){
    s_sample=strsplit(sample,"_")
    patient=s_sample[[1]][1]
    sort=s_sample[[1]][2]
    patients.big=c(patients.big,patient)
    sorts.big=c(sorts.big,sort)
  }
  scrna.big$patient=patients.big
  scrna.big$sort=sorts.big
  return(scrna.big)
}

# Merge other groups / indiv dataset attributes - GENERAL USE
add_groups<-function(scrna.big,ngroup_levels,has_splitsamp,replace_origident){
  samples.big=scrna.big$orig.ident
  samples_compiled.big=c()
  patients.big=c()
  sorts.big=c()
  if(ngroup_levels==3){
    patient_timepoints=c()
    timepoints=c()
  }
  for(sample in samples.big){
    s_sample=strsplit(sample,"_")
    patient=s_sample[[1]][1]
    sort=s_sample[[1]][2]
    patients.big=c(patients.big,patient)
    sorts.big=c(sorts.big,sort)
    if(ngroup_levels==3){
      timepoint=s_sample[[1]][3]
      timepoints=c(timepoints,timepoint)
      if(has_splitsamp==TRUE){
        splitsamp<-s_sample[[1]][4]# this is if the slice of data has multiple  files
      }
      sample_compiled<-paste(patient,sort,timepoint,sep="_")
      patient_timepoint<-paste(patient,timepoint,sep="_")
      patient_timepoints<-c(patient_timepoints,patient_timepoint)
    }
    else{
      if(has_splitsamp==TRUE){
        splitsamp<-s_sample[[1]][3]# this is if the slice of data has multiple  files
      }
      sample_compiled<-paste(patient,sort,sep="_")
    }
    samples_compiled.big=c(samples_compiled.big,sample_compiled)
    
  }
  scrna.big$patient=patients.big
  scrna.big$sort=sorts.big
  scrna.big$samples_compiled=samples_compiled.big
  if(replace_origident==TRUE){
    scrna.big$origident=samples_compiled.big
  }
  if(ngroup_levels==3){
    scrna.big$patient_timepoint<-patient_timepoints
    scrna.big$timepoint<-timepoints
  }
  return(scrna.big)
}



# Grab slice of seurat object matching groupby=groupname condition
get_groups<-function(scrna,groupby,groupname){
  if(groupby=='sort'){
    scrna.pt2<-scrna[,scrna$sort==groupname]
  }
  if(groupby=='sample'){
    scrna.pt2<-scrna[,scrna$sample==groupname]
  }
  if(groupby=='patient'){
    scrna.pt2<-scrna[,scrna$patient==groupname]
  }
  if(groupby=='timepoint'){
    scrna.pt2<-scrna[,scrna$timepoint==groupname]
  }
  return(scrna.pt2)
}


# Returns a vector with all cells that have > 0 (background) expression of a given marker gene
# Need to modify this post-20200714 to  return different levels of exp
get_nonzero_marker_cells<-function(gex,marker,pos){#gex is GetAssayData(scrna)
  gex_marker=gex[rownames(gex)==marker,]
  if(length(gex_marker)==0){
    return(c(paste(marker,"negative",sep="-")))
  }
  lower_bound <- mean(gex_marker) + sd(gex_marker)# define boundary for empty cells (cells with <  1SD ABOVE MEAN)
  # Filter out empty cells (if pos==TRUE) ; filter out non-empty cells (if pos==FALSE)
  if(pos){
    gex_marker=gex_marker[gex_marker>lower_bound]
  }
  else{
    gex_marker=gex_marker[gex_marker<=lower_bound]
  }
  return(gex_marker)
}

# Takes a vector of cellids and return a seurat-object-wide membership vector
# consisting of labels for in- vs out-class membership
cellids_membership<-function(scrna,cellids,inclass,outclass){
  all_cellids<-colnames(scrna)
  membership<-all_cellids%in%cellids
  membership<-ifelse(membership==TRUE,inclass,outclass)
  return(membership)
}


# Takes multiple membership vectors and combines them, replacing nacol cells
# and merging where (non-nacol) overlap occurs
merge_cellids_membership<-function(memberships,nacol){
  l=length(memberships[[1]])
  merged_membership=c()
  for(i in 1:l){
    haslabels=c()
    for(membership in memberships){
      if(!membership[i]==nacol){
        haslabels<-c(haslabels,membership[i])
      }
    }
    if(length(haslabels)>0){
      haslabels<-paste(haslabels,collapse="+")
    }
    else{
      haslabels=nacol
    }
    merged_membership<-c(merged_membership,haslabels)
  }
  return(merged_membership)
}