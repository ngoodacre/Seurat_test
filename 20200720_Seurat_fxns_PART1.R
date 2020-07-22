
####################################################################################################
####################################################################################################
####################################################################################################
####################################          PART I                     ###########################
####################################    Ingests h5 files                 ###########################
####################################    collapses by patient/splitby     ###########################
####################################    normalizes and saves             ###########################
####################################################################################################
####################################################################################################
####################################################################################################


####################################################################################################
####################################################################################################
#### COMMAND BLOCK FOR PART 1
RUN_PART_1<-function(datadir,bad_files,refgenes,ngroup_levels,splitby,session_name1,figs_folder,date){
  message(paste("Retrieving h5 files from",datadir,sep=" "))
  h5files=list.files(datadir)[5:8]
  print(getwd())
  message("Setting up Seurat objects")
  r<-setup_seurat_objects(datadir,h5files,bad_files,ngroup_levels,analysis)
  seurat_objects<-r[[1]]
  samples<-r[[2]]
  
  #!message("Patching in any missing referene gene info (as zero vectors)")
  #!seurat_objects<-patch_gene_rows(seurat_objects,refgenes)
  message("Merging Seurat objects")
  scrna.merged<-merge_seurat_objects(seurat_objects,samples,analysis,ngroup_levels)
  
  message("Splitting merged Seurat object by",splitby,sep=" ")
  scrna.bysplit<-split_merged_seurat_object(scrna.merged,splitby,figs_folder,date)
  f=1
  splits<-names(scrna.bysplit)
  for(scrna.split in scrna.bysplit){
    split<-splits[f]
    session_name1_split<-paste(sub(".RDS","",session_name1),paste("split=",split,".RDS",sep=""),sep="_")
    message(paste("Saving PART1 session split as",session_name1_split,sep=" "))
    saveRDS(scrna.split, file = session_name1_split)
    f=f+1
  }
  #return(scrna.bysplit)
}


####################################################################################################
####################################################################################################
#### SET UP THE SEURAT OBJECT
# Load the raw (non-normalized data).
#scrna.data <- Read10X_h5(data.dir = getwd())
setup_seurat_objects<-function(datadir,h5files,badfiles,ngroups,analysis){
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
    scrna.data <- Read10X_h5(paste(datadir,h5file,sep="/"))
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
####################################################################################################
#### UPDATE ALL SEURAT OBJECTS IN LIST TO CONTAIN REF/PROVIDED LIST (WITH ZEROS)
patch_gene_rows<-function(seurat_object_list,refgenes){
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
#### MERGED SEURAT OBJECTS (USUALLY PRE-INTEGRATION)
merge_seurat_objects1<-function(seurat_objects,samples,analysis,ngroup_levels){
  scrna.merged=merge(seurat_objects[[1]], y =c(seurat_objects[[2]],seurat_objects[[3]],seurat_objects[[4]],seurat_objects[[5]],
                                               seurat_objects[[6]],seurat_objects[[7]],seurat_objects[[8]],seurat_objects[[9]],seurat_objects[[10]]),
                                               add.cell.ids = samples, project = "scrna")
  #### Merge other groups / indiv dataset attributes
  scrna.merged<-add_groups(scrna.merged,ngroup_levels,FALSE,FALSE)
  return(scrna.merged)
}
merge_seurat_objects<-function(seurat_objects,samples,analysis,ngroup_levels){
  scrna.merged=merge(seurat_objects[[1]], y =c(seurat_objects[[2]],seurat_objects[[3]],seurat_objects[[4]]),
                     add.cell.ids = samples, project = "scrna")
  #### Merge other groups / indiv dataset attributes
  scrna.merged<-add_groups(scrna.merged,ngroup_levels,FALSE,FALSE)
  return(scrna.merged)
}


####################################################################################################
####################################################################################################
#### RE-SPLIT MERGED DATA - by splitby param (default='patient') 
split_merged_seurat_object<-function(scrna.merged,splitby,figs_folder,date){
  scrna.bysplit<-SplitObject(scrna.merged,split.by=splitby)
  f=1
  for(scrna.split in scrna.bysplit){
    split=names(scrna.bysplit)[f]
    scrna.bysplit[split]<-QC_norm(scrna.split,figs_folder,split,date)
    f=f+1
  }
  return(scrna.bysplit)
}


####################################################################################################
#### QUALITY CELL/GENE CONTROL AND NORMALIZATION
QC_norm <- function(scrna,figs_folder,name,date){
  print("Performing QC and normalization")
  scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^MT-")
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
  ### implement filtering to remove empty cells, multiplet cells (heuristic from Monocle Vignette at: http://cole-trapnell-lab.github.io/monocle-release/docs/)
  scrna.ncounts=scrna$nCount_RNA
  upper_bound <- 10^(mean(log10(scrna.ncounts)) + 2*sd(log10(scrna.ncounts)))# define multiplets
  lower_bound <- 10^(mean(log10(scrna.ncounts)) - 2*sd(log10(scrna.ncounts)))# define empty cells
  scrna=scrna[,scrna.ncounts>lower_bound&scrna.ncounts<upper_bound]
  #### NORMALIZING THE DATA
  scrna <- NormalizeData(scrna, normalization.method = "LogNormalize", scale.factor = 10000) #specify defaults
  scrna <- PercentageFeatureSet(scrna, pattern = "^MT-", col.name = "percent.mt")
  return(scrna)
}