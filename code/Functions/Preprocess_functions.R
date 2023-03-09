# (1) object to sce ----
# Devuelve un objeto SCE listo para la funciones ldf y sce.markers 
# Inputs:
#  xx: puede ser un SCE, Seurat o una matriz de counts
#  assay.name.raw: el nombre del assay con la matriz de raw counts
#  assay.name.normalization: el nombre del assay con la matriz normalizada a usar en el shiny,
#  metadata: un dataframe con la metadata de las celulas
#  toFactors: el nombre de las columnas de la metadata o colData para pasar a factores
#  chosen.hvg: lista de nombres de los genes mas variables
#  nHVGs: si chosen.hvg=NULL calculo los nHVGs genes mas variables
#  nPCs: número de componentes principales a calcular para PCA
#  verbose: que la función de información de los pases que está haciendo
#  calcRedDim: si tiene que computar las dimensiones reducidas (PCA, UMAP, TSNE, UMAP2D, TSNE2D)
createSCEobject <- function(xx,
                            assay.name.raw="counts",
                            assay.name.normalization="logcounts",
                            metadata=NULL,
                            toFactors=NULL,
                            chosen.hvg=NULL,
                            nHVGs=3000,
                            nPCs=50,
                            calcRedDim=TRUE,
                            verbose=TRUE){
  
  csceo <- list()
  
  library(scater)
  library(scran)
  library(Matrix)
  library(SingleCellExperiment)
  library(doParallel)
  
  ##Check if there are colnames and rownames in the objects 
  if(is.null(rownames(xx)) | is.null(colnames(xx))){
    stop('Missing row names or column names.')
  }
  #Check if the user put twice the same partition.
  toFactors <- unique(toFactors)
  
  if(is.null(toFactors)){
    warning('No toFactors specified, a quick clusterization will be computed.')
    toFactors <- "scx.clust"
  }
  #If isn't a Seurat or a SCE object, and the metadata is null (because the metadata could be inside the object)
  #The toFactors become the default.
  if(is.null(metadata) & class(xx)[1]!="Seurat" & class(xx)[1]!="SingleCellExperiment"){
    warning('No metadata specified, a quick clusterization will be computed.')
    toFactors <- "scx.clust"
  }
  if(verbose) cat('Creating SCE object...')
  # Seurat to sce ----
  if(class(xx)[1]=="Seurat"){
    # changing assay.name parameters because as.SingleCellExperiment fills 'counts' and 'logcounts' assays
    warning("Seurat object is detected, 'assay.name.raw' & 'assay.name.normalization' will be set to default.")
    assay.name.raw <- "counts"
    assay.name.normalization <- "logcounts"
    
    library(Seurat)
    xx.sce <- as.SingleCellExperiment(xx) 
    if((all(toFactors!="scx.clust")) & (!all(toFactors %in% names(colData(xx.sce))))){
      warning('at least one factor is not present in metadata')
    }
    if(nrow(xx@assays$RNA@meta.features)>0){
      rowData(xx.sce) <- xx@assays$RNA@meta.features
    }
    if(length(VariableFeatures(xx))>0){
      chosen.hvg <- VariableFeatures(xx)
    }
    # matrix + metadata to sce ----
  } else if (class(xx)[1] %in% c("dgCMatrix", "Matrix", "matrix")){
    xx.sce <- SingleCellExperiment(list(counts=xx))
    if(!is.null(metadata)){
      colData(xx.sce) <- c(colData(xx.sce), metadata)
      if((!("scx.clust" %in% toFactors)) & (!all(toFactors %in% names(colData(xx.sce))))){
        warning('at least one factor is not present in metadata')
      }
    } 
    # else {
    #   stop('metadata not found')
    # }
  } else if (class(xx)[1]=="SingleCellExperiment"){
    xx.sce <- xx
    if(!is.null(metadata)){
      colData(xx.sce) <- c(colData(xx.sce),metadata)
      if((toFactors!="scx.clust") & (!all(toFactors %in% names(colData(xx.sce))))){
        warning('at least one factor is not present in metadata')
      }
    } 
    # if(!sum(toFactors %in% names(colData(xx.sce)))==length(toFactors)){
    #   stop('at least one factor is not present in metadata')
    # }
  } else {
    stop('xx must be an object of the class Seurat, SingleCellExperiment or Matrix')
  }
  if(verbose) cat(' Finished\n')
  
  #If there are not an assay called counts, get the first assay from the list and rename it as counts (to be consistent after on the code)
  # if(!('counts' %in% names(assays(xx.sce)))){
  #   names(assays(xx.sce))[1] <- 'counts'
  # }
  
  # QC ----
  if(verbose) cat('Computing QC metrics...')
  if(!assay.name.raw %in% names(assays(xx.sce))){
    stop(paste0('Assay ',name.assay.raw,' not found in SCE object'))
  }
  nCounts <- apply(assay(xx.sce, assay.name.raw),2,sum)
  nFeatures <- apply(assay(xx.sce, assay.name.raw),2,function(x){sum(x>0)})
  df.qc <- data.frame(row.names = colnames(xx.sce), nCounts = nCounts, nFeatures = nFeatures)
  colData(xx.sce) <- c(colData(xx.sce),df.qc)
  if(verbose) cat(' Finished\n')
  
  
  # Normalization ----
  clust <- quickCluster(xx.sce, assay.type = assay.name.raw)
  xx.sce$scx.clust <- clust
  if(!assay.name.normalization %in% names(assays(xx.sce))){
    if(verbose) cat('Computing normalization...')
    set.seed(123457)
    clust <- quickCluster(xx.sce, assay.type = assay.name.raw)
    xx.sce$scx.clust <- clust
    xx.sce <- computeSumFactors(xx.sce,cluster=clust,min.mean=0.1, assay.type = assay.name.raw)
    # xx.sce <- scuttle::pooledSizeFactors(xx.sce,cluster=clust,min.mean=0.1, assay.type = assay.name.raw)
    xx.sce <- logNormCounts(xx.sce, assay.type = assay.name.raw, name="logcounts")
    if(verbose) cat(' Finished\n')
  }
  
  
  # HVGs ----
  if(is.null(chosen.hvg)){
    if(verbose) cat('Computing HVGs...')
    mgv <- modelGeneVar(xx.sce,span=.8, assay.type = assay.name.normalization)
    rowData(xx.sce) <- cbind(rowData(xx.sce), hvg.mvBio=mgv$bio)
    chosen.hvg <- rank(-rowData(xx.sce)$hvg.mvBio) <= 3000 & rowData(xx.sce)$hvg.mvBio>0
    chosen.hvg <- rownames(xx.sce)[chosen.hvg]
    if(verbose) cat(' Finished\n')
  }
  # chosen.hvg <- chosen.hvg[!chosen.hvg%in%c("Ehd2","Espl1","Jarid1d","Pnpla4","Rps4y1","Xist","Tsix",
  #                "Eif2s3y", "Ddx3y", "Uty","Kdm5d","Rpl26","Gstp1","Rpl35a",
  #                "Erh","Slc25a5","Pgk1","Eno1","Tubb2a","Emc4", "Scg5")]
  
  
  # red dims ----
  #If the user specify no calcRedDim but:
  #there is no dimRed calculated, it calcules all
  #there is no 2D or 3D redDim caclulated, it calcules all
  #there is no 2D, it calculates only the 2D.
  #there is no 3D, it calculates only the 3D.
  
  runDim <- c("PCA", "TSNE", "UMAP", "TSNE2D", "UMAP2D")
  if(!calcRedDim){
    if(length(reducedDimNames(xx.sce))<1 | all(!(sapply(reducedDims(xx.sce),ncol) %in% c(2,3)))){
      calcRedDim <- TRUE
      runDim <- c("PCA", "TSNE", "UMAP", "TSNE2D", "UMAP2D")
    } else if(all(sapply(reducedDims(xx.sce),ncol) != 2)){
      calcRedDim <- TRUE
      runDim <- c("PCA","TSNE2D", "UMAP2D")
    } else if(all(sapply(reducedDims(xx.sce),ncol) != 3)){
      calcRedDim <- TRUE
      runDim <- c("PCA", "TSNE", "UMAP")
    }
  }
  
  if(calcRedDim){
    if(verbose) cat('Computing following reduced dims:',paste(runDim, collapse = ' '),'\n')
    xx.sce <- applyReducedDim(xx.sce, runDim, chosen.hvg, nPCs, assay.name.normalization, prefix.name="SCX_", verbose)
    if(verbose) cat('Finished\n')
  }
  
  # changing name for normalizated matrix to be inputed in shiny app
  if(assay.name.normalization!="logcounts" & "logcounts" %in% assayNames(xx.sce)){
    assay(xx.sce, "logcounts") <- NULL
  }
  assayNames(xx.sce)[which(assayNames(xx.sce)==assay.name.normalization)] <- "logcounts"
  
  if(verbose) cat('Computing logcounts normalized...')
  #logcounts normalized ----
  if(!"logcounts.norm" %in% names(assays(xx.sce))){
    assays(xx.sce)$logcounts.norm <- t(apply(assay(xx.sce, "logcounts"), 1, function(x){x/max(x)}))
    assays(xx.sce)$logcounts.norm <- Matrix(assays(xx.sce)$logcounts.norm, sparse = TRUE)
  }
  if(verbose) cat(' Finished\n')
  
  
  #to factors ----
  if(verbose) cat('Changing factors from toFactors...')
  #Check the factors that are in the colData.
  tfs <- setNames(nm=toFactors,object = toFactors %in% names(colData(xx.sce)))
  #Check if this columns have more than one level.
  if(sum(tfs) > 0){
    vld <- sapply(colData(xx.sce)[,toFactors[tfs],drop=F],function(x){length(unique(x))>1})   
    tfs[names(vld)] <-  tfs[names(vld)] & vld
  }
  
  ttoFactors <- names(tfs)[tfs]
  if(sum(!tfs)>0) warning("Can't find ",paste0(names(tfs)[!tfs],collapse = ' & ')," in metadata or they have only one level")
  if(length(ttoFactors) == 0){
    warning("Any factor passed the controls, add the quick cluster as toFactors")
    ttoFactors <- "scx.clust"
  }
  
  #Transform all the character columns to factor to be able to be selected in the shinyApp. 
  cols <- sapply(colData(xx.sce), function(x){(is.character(x) | is.factor(x)) & length(unique(x)) < 30})
  if(any(cols)){
    colData(xx.sce)[,cols] <- lapply(colData(xx.sce)[,cols,drop=F], function(x){
        x <- as.factor(x)
        x <- droplevels(x)
        x
      })
  if(verbose) cat(' Finished\n')
  
  
  # Attaching SCE to output ----
  csceo[["SCE"]] <- xx.sce
  
  
  if(verbose) cat('Computing sce.markers...')
  # sce markers ----
  sce.markers <- list()
  for(i in ttoFactors){
    sce.markers[[i]] <- findMarkers(xx.sce, 
                                    assay.type = "logcounts",
                                    group = colData(xx.sce)[,i],
                                    direction="any",pval.type="all",log.p=T,full.stats=T)
  }
  if(verbose) cat(' Finished\n')

  # Attaching sce.markers to output
  csceo[["sce.markers"]] <- sce.markers
  
  
  if(verbose) cat('Computing ldf function\n')
  # ldf functions ----
  ldf <- list()
  for(i in ttoFactors){
    ldf[[i]] <- ldf_func(xx.sce, i)
  }
  if(verbose) cat('Finished\n')
  
  
  # Attaching ldf to output ----
  csceo[["ldf"]] <- ldf
  
  gc() # free memory
  return(csceo)
}

# (2) Red Dim ----
applyReducedDim <- function(sce, reddimstocalculate, chosen.hvgs, nPCs, assayname, prefix.name="SCX_",verbose=TRUE){
  namepca <- paste0(prefix.name,"PCA")
  if("PCA"%in%reddimstocalculate){
    if(verbose) cat("\t Calculating PCA...")
    set.seed(12534)
    sce <- runPCA(sce, subset_row=chosen.hvgs, ncomponents=nPCs, name=namepca, exprs_values=assayname)
    if(verbose) cat(' Finished','\n')
  }
  if("TSNE"%in%reddimstocalculate){
    if(verbose) cat("\t Calculating TSNE...")
    set.seed(1111011)
    sce <- runTSNE(sce,dimred=namepca,n_dimred=20,ncomponents=3,name=paste0(prefix.name,"TSNE"),exprs_values=assayname)
    if(verbose) cat(' Finished','\n')
  }
  if("UMAP"%in%reddimstocalculate){
    if(verbose) cat("\t Calculating UMAP...")
    set.seed(1111011)
    sce <- runUMAP(sce,dimred=namepca,n_dimred=20,ncomponents=3,name=paste0(prefix.name,"UMAP"),exprs_values=assayname)
    if(verbose) cat(' Finished','\n')
  }
  if("TSNE2D"%in%reddimstocalculate){
    if(verbose) cat("\t Calculating TSNE2D...")
    set.seed(1111011)
    sce <- runTSNE(sce,dimred=namepca,name=paste0(prefix.name,"TSNE2D"),exprs_values=assayname)
    if(verbose) cat(' Finished','\n')
  }
  if("UMAP2D"%in%reddimstocalculate){
    if(verbose) cat("\t Calculating UMAP2D...")
    set.seed(1111011)
    sce <- runUMAP(sce,dimred=namepca,name=paste0(prefix.name,"UMAP2D"),exprs_values=assayname)
    if(verbose) cat(' Finished','\n')
  }
  
  return(sce)
}

# (3) ldf func ----
ldf_func <- function(sce, partition,minSize=50){
  
  cat(partition, ":\n", sep = "")
  # # le saco los factores vacios
  # colData(sce)[,partition] <- droplevels(colData(sce)[,partition])
  
  # calculo lfmrk ----
  type           <- c('all','any','some')
  lfmrk <- list()
  for(itype in seq_along(type)){
    lfmrk[[type[itype]]]    <- findMarkers(sce,
                                           assay.type = "logcounts",
                                           groups=colData(sce)[,partition],
                                           full.stats=TRUE,
                                           direction='up',
                                           pval.type=type[itype])
  }
  
  
  # calculo ldf ----
  ldf_t <-list()
  numCores <- detectCores() - 1
  registerDoParallel(numCores)
  # lfmrk <- lfmrk_t
  
  lab   <- colData(sce)[,partition]
  tt   <- table(lab)
  labOK <- tt > minSize
  labOK <- names(labOK)[labOK]
  
  lfmrk <- lapply(lfmrk,function(x){x[labOK]})
  
  if(length(lfmrk[[1]])>0){ #At least one group to calculate
    for(ic in seq_along(lfmrk[[1]])){
      coi <- names(lfmrk[[1]])[ic]
      cat(coi,'\n')
      any  <- rownames(lfmrk[['any']][[coi]])[lfmrk[['any']][[coi]][,'FDR']<0.05 & 
                                                lfmrk[['any']][[coi]][,'Top']<=10]
      all  <-rownames(lfmrk[['all']][[coi]])[lfmrk[['all']][[coi]][,'FDR']<0.05]
      some <-rownames(lfmrk[['some']][[coi]])[lfmrk[['some']][[coi]][,'FDR']<0.05]
      u <- unique(c(any,all,some))
      maux <- cbind(all=u%in%all,any=u%in%any,some=u%in%some)
      rownames(maux)<-u
      a     <- apply(maux,1,function(x){sum(x*c(4,2,1))})
      mrkrs <- sort(a,decreasing=TRUE)
      a4    <- names(mrkrs)[mrkrs>=1]
      
      # (6.1.1) Calculo de Boxcor ----
      
      cat("\t Calculating boxcor\n")
      
      if(length(a4) > 0){
        Z   <- assay(sce, "logcounts")[a4,]  
        
        pattern <- rep(0,ncol(Z))
        pattern[colData(sce)[,partition]%in%coi] <- 1
        boxcor  <-apply(Z,1,
                        function(x){
                          return(cor(x,pattern))
                        }
        )
        df <- data.frame(boxcor=boxcor,
                         robustness=a[names(mrkrs)]
        )
      }
      else {
        df <- NULL
      }
      ldf_t[[coi]]<-df
    }
  }
  return(ldf_t)
}
