# Object to sce ----
#' SingleCellExperiment object ready for use in scX
#'
#' `createSCEobject()` returns a list which includes a SingleCellExperiment object with counts, at
#'		least one clusterization, gene markers for each clusterization, reduced dimensions for
#'		visualization, and any additional data provided. This list is used as input when launching
#'		the scX app.
#' 
#' @param xx Either a matrix with counts, or a SCE or Seurat object.
#' @param assay.name.raw Assay name for raw counts matrix if object is a SCE. Defaults to `counts`.
#' @param assay.name.normalization Assay name for normalized matrix if present in SCE. If not present,
#' 		it computes `logcounts` (default).
#' @param metadata Optional dataframe with cell metadata.
#' @param partitionVars Optional metadata column names (or `colData`) to use as factors. If `NULL` (default),
#'		a quick clusterization will be computed.
#' @param metadataVars Additional metadata column names (or `colData`) to use only for coloring in plots.
#'    If `NULL` (default), only `partitionVars` columns will be available for coloring plots.
#' @param chosen.hvg Optional list of Highly Variable Genes.
#' @param nHVGs Number of Highly Variable Genes to use if `chosen.hvg=NULL`. Defaults to 3000.
#' @param nPCs Number of Principal Components to use in PCA. Defaults to 50.
#' @param verbose Logical for step by step status while function is running. Defaults to `TRUE`.
#' @param calcRedDim Logical indicating whether to compute reduced dimensions (PCA, UMAP, TSNE, UMAP2D, 
#'		TSNE2D). Defaults to `TRUE`.
#' @param paramFindMarkers List of parameters to pass to `scran::findMarkers(...)` to compute marker 
#'		genes for clusters. Defaults to `list(test.type="wilcox", pval.type="all", direction="any")`.
#' @param calcAllPartitions Logical indicating whether to force the computation of markers and DEGs from
#'		the entire list of `partitionVars`.
#' @param cells2keep List of names for cells to keep when subsampling data. NOTE: Subsampling is only 
#'		used in large datasets for visual purposes, and it does not affect computations.
#' @param descriptionText Optional short description of the object being analized. This can help when 
#'		working with multiple tabs. 
#' @returns List with a SingleCellExperiment object and additional data ready for use in `scX`.
#' @export
createSCEobject <- function(xx,
                            assay.name.raw="counts",
                            assay.name.normalization="logcounts",
                            metadata=NULL,
                            partitionVars=NULL,
                            metadataVars=NULL,
                            chosen.hvg=NULL,
                            nHVGs=3000,
                            nPCs=50,
                            calcRedDim=TRUE,
                            paramFindMarkers=list(test.type="wilcox", pval.type="all", direction="any"),
                            calcAllPartitions=FALSE,
                            cells2keep=NULL,
                            descriptionText=NULL,
                            verbose=TRUE){
  
  csceo <- list()
  ##Check if there are colnames and rownames in the objects 
  if(is.null(rownames(xx)) | is.null(colnames(xx))){
    stop('Missing row names or column names.')
  }
  #Check for repeated partitions
  partitionVars <- unique(partitionVars)
  
  if(is.null(partitionVars)){
    warning('No partitionVars specified, a quick clusterization will be computed.')
    partitionVars <- "scx.clust"
  }
  #If it isn't a Seurat or SCE object, and metadata is null (metadata could be inside the object)
  #The partitionVars become the default.
  if(is.null(metadata) & class(xx)[1]!="Seurat" & class(xx)[1]!="SingleCellExperiment"){
    warning('No metadata specified, a quick clusterization will be computed.')
    partitionVars <- "scx.clust"
  }
  if(verbose) cat('Creating SCE object...')
  # Seurat to sce ----
  if(class(xx)[1]=="Seurat"){
    # changing assay.name parameters because as.SingleCellExperiment fills 'counts' and 'logcounts' assays
    warning("Seurat object is detected, 'assay.name.raw' & 'assay.name.normalization' will be set to default.")
    assay.name.raw <- "counts"
    assay.name.normalization <- "logcounts"
    
    xx.sce <- Seurat::as.SingleCellExperiment(xx) 
    if((all(partitionVars!="scx.clust")) & (!all(partitionVars %in% names(colData(xx.sce))))){
      warning('at least one partition is not present in metadata')
    }
    if(nrow(xx@assays$RNA@meta.features)>0){
      rowData(xx.sce) <- xx@assays$RNA@meta.features
    }
    if(length(Seurat::VariableFeatures(xx))>0){
      chosen.hvg <- Seurat::VariableFeatures(xx)
    }
    # matrix + metadata to sce ----
  } else if (class(xx)[1] %in% c("dgCMatrix", "Matrix", "matrix")){
    if(is.null(rownames(xx)) | is.null(colnames(xx))){
      stop('Matrix must have rownames "gene_id" and colnames "barcodes"')
    }
    xx.sce <- SingleCellExperiment(list(counts=xx))
    if(!is.null(metadata)){
      if(all(colnames(xx.sce) %in% rownames(metadata))){
        colData(xx.sce) <- cbind(colData(xx.sce), metadata[colnames(xx.sce),,drop=F])
        if((!("scx.clust" %in% partitionVars)) & (!all(partitionVars %in% names(colData(xx.sce))))){
          warning('at least one partition is not present in metadata')
        }
      } else {
        stop('Some cells in metadata are not present in colnames matrix.\n')
        }
    }
    # else {
    #   stop('metadata not found')
    # }
  } else if (class(xx)[1]=="SingleCellExperiment"){
    xx.sce <- xx
    if(!is.null(metadata)){
      colData(xx.sce) <- cbind(colData(xx.sce),metadata)
      if((partitionVars!="scx.clust") & (!all(partitionVars %in% names(colData(xx.sce))))){
        warning('at least one partition is not present in metadata')
      }
    } 
    # if(!sum(partitionVars %in% names(colData(xx.sce)))==length(partitionVars)){
    #   stop('at least one factor is not present in metadata')
    # }
  } else {
    stop('xx must be an object of the class Seurat, SingleCellExperiment or Matrix')
  }
  if(verbose) cat(' Finished\n')
  
  #to factors ----
  if(verbose) cat('Changing factors from partitionVars...')
  #Check the factors that are in the colData.
  tfs <- setNames(nm=partitionVars,object = partitionVars %in% names(colData(xx.sce)))
  #Check if columns have more than one level.
  if(sum(tfs) > 0){
    vld <- sapply(colData(xx.sce)[,partitionVars[tfs],drop=F],function(x){length(unique(x))>1})   
    tfs[names(vld)] <-  tfs[names(vld)] & vld
  }
  
  ttoFactors <- names(tfs)[tfs]
  if(sum(!tfs)>0) warning("Can't find ",paste0(names(tfs)[!tfs],collapse = ' & ')," in metadata or they have only one level")
  if(length(ttoFactors) == 0){
    warning("No partition passed the controls, a quick clusterization will be computed.")
    ttoFactors <- "scx.clust"
  }else{
	colData(xx.sce)[ttoFactors] <- lapply(colData(xx.sce)[ttoFactors], as.factor) # in case a partition is numeric
  }
  
  #Transform all the character columns to factor to be able to be selected in the shinyApp. 
  cols <- sapply(colData(xx.sce), function(x){(is.character(x) | is.factor(x)) & length(unique(x)) <= 30})
  if(any(cols)){
    colData(xx.sce)[,cols] <- lapply(colData(xx.sce)[,cols,drop=F], function(x){
        x <- as.factor(x)
        x <- droplevels(x)
        x
      })
	}
  if(verbose) cat(' Finished\n')
    
  # Checking number of levels of ttoFactors for calculations ----
  if(!calcAllPartitions){
    allToFactors <- sapply(ttoFactors, function(x){length(unique(colData(xx.sce)[,x]))})>30
    if(all(allToFactors)){
      stop(paste0(paste0(names(allToFactors)[allToFactors], collapse =  ' & ')," has more than 30 levels. If you want to compute it anyway set 'calcAllPartitions' as TRUE"))
    } else if(any(allToFactors)) {
      warning(paste0(paste0(names(allToFactors)[allToFactors], collapse =  ' & ')," has more than 30 levels. They wont be used to compute markers and DEGs. If you want to compute it anyway set 'calcAllPartitions' as TRUE"))
      ttoFactors <- ttoFactors[!allToFactors]
    }
  }
  
  #If there are not an assay called counts, get the first assay from the list and rename it as counts (to be consistent after on the code)
  # if(!('counts' %in% names(assays(xx.sce)))){
  #   names(assays(xx.sce))[1] <- 'counts'
  # }
  
  # QC ----
  if(verbose) cat('Computing QC metrics...')
  if(!assay.name.raw %in% names(assays(xx.sce))){
    if(assay.name.normalization %in% names(assays(xx.sce))){
      warning(paste0('Assays ',assay.name.raw,' not found in SCE object'))
      xx.sce$nCounts   <- NA
      xx.sce$nFeatures <- NA
    } else {
      stop(paste0('Assay ',paste(assay.name.raw, assay.name.normalization, sep = ' & '),' not found in SCE object'))
    }
  } else {
      xx.sce$nCounts <- colSums(assay(xx.sce, assay.name.raw))
      if(class(assay(xx.sce, assay.name.raw))[1]%in%c("dgCMatrix")){
          xx.sce$nFeatures <- diff(assay(xx.sce, assay.name.raw)@p)#apply(assay(xx.sce, assay.name.raw),2,function(x){sum(x>0)})
      } else {
          xx.sce$nFeatures <- apply(assay(xx.sce, assay.name.raw),2,function(x){sum(x>0)})
      }
  }
  if(verbose) cat(' Finished\n')
  
  
  # Normalization ----
  if(!assay.name.normalization %in% names(assays(xx.sce))){
    if(assay.name.raw %in% names(assays(xx.sce))){
      if(verbose) cat('Computing normalization...')
      set.seed(123457)
      clust <- scran::quickCluster(xx.sce, assay.type = assay.name.raw)
      xx.sce$scx.clust <- clust
      xx.sce <- scran::computeSumFactors(xx.sce,cluster=clust,min.mean=0.1, assay.type = assay.name.raw)
      xx.sce <- scater::logNormCounts(xx.sce, assay.type = assay.name.raw, name="logcounts")
      if(verbose) cat(' Finished\n')
    }
    # xx.sce <- scuttle::pooledSizeFactors(xx.sce,cluster=clust,min.mean=0.1, assay.type = assay.name.raw)
  } else if ( "scx.clust" %in% partitionVars ) {
      if(verbose) cat('Computing clusters...')
      clust <- scran::quickCluster(xx.sce, assay.type = assay.name.normalization)
      xx.sce$scx.clust <- clust
      if(verbose) cat(' Finished\n')
  }
  
  
  # HVGs ----
  if(is.null(chosen.hvg)){
    if(verbose) cat('Computing HVGs...')
    mgv <- scran::modelGeneVar(xx.sce,span=.8, assay.type = assay.name.normalization)
    rowData(xx.sce) <- cbind(rowData(xx.sce), hvg.mvBio=mgv$bio)
    chosen.hvg <- rank(-rowData(xx.sce)$hvg.mvBio) <= 3000 & rowData(xx.sce)$hvg.mvBio>0
    chosen.hvg <- rownames(xx.sce)[chosen.hvg]
    if(verbose) cat(' Finished\n')
  }
  # chosen.hvg <- chosen.hvg[!chosen.hvg%in%c("Ehd2","Espl1","Jarid1d","Pnpla4","Rps4y1","Xist","Tsix",
  #                "Eif2s3y", "Ddx3y", "Uty","Kdm5d","Rpl26","Gstp1","Rpl35a",
  #                "Erh","Slc25a5","Pgk1","Eno1","Tubb2a","Emc4", "Scg5")]
  
  
  # Reduced dimensions ----
  #If user specifies no calcRedDim but:
  #there is no dimRed calculated -> it calculates all
  #there is no 2D or 3D redDim caclulated -> it calculates all
  #there is no 2D -> it calculates only 2D
  #there is no 3D -> it calculates only 3D
  
  runDim <- c("PCA", "TSNE", "UMAP", "TSNE2D", "UMAP2D")
  if(!calcRedDim){
    if(length(reducedDimNames(xx.sce))<1 | all(!(sapply(reducedDims(xx.sce),ncol) %in% c(2,3)))){
      calcRedDim <- TRUE
      #runDim <- c("PCA", "TSNE", "UMAP", "TSNE2D", "UMAP2D")
    } else if(all(sapply(reducedDims(xx.sce),ncol) != 2)){
      calcRedDim <- TRUE
      runDim <- c("PCA","TSNE2D", "UMAP2D")
    } else if(all(sapply(reducedDims(xx.sce),ncol) != 3)){
      calcRedDim <- TRUE
      runDim <- c("PCA", "TSNE", "UMAP")
    }
  }
  
  if(calcRedDim){
    if(verbose) cat('Computing the following reduced dims:',paste(runDim, collapse = ' '),'\n')
    xx.sce <- applyReducedDim(xx.sce, runDim, chosen.hvg, nPCs, assay.name.normalization, prefix.name="SCX_", verbose)
    if(verbose) cat('Finished\n')
  }
  
  # changing name for normalized matrix to be input in shiny app
  if(assay.name.normalization!="logcounts" & "logcounts" %in% assayNames(xx.sce)){
    assay(xx.sce, "logcounts") <- NULL
  }
  assayNames(xx.sce)[which(assayNames(xx.sce)==assay.name.normalization)] <- "logcounts"
  
  if(verbose) cat('Computing logcounts normalized...')
  #logcounts normalized ----
  if(!"logcounts.norm" %in% names(assays(xx.sce))){
		sparse_mat <- as(assay(xx.sce, "logcounts"), "sparseMatrix")
		row_maxs <- qlcMatrix::rowMax(sparse_mat)
		maxdiag <- Diagonal(x = 1/as.vector(row_maxs))
		scaled_sparse <- maxdiag %*% sparse_mat
		rownames(scaled_sparse) <- rownames(sparse_mat)

		assays(xx.sce)$logcounts.norm <- scaled_sparse
  }
  if(verbose) cat(' Finished\n')
  

  # Subsetting SCE object ----
  # Keep only colData names that will be use for coloring plots
  if(!is.null(metadataVars)){
    coldatanames <- names(colData(xx.sce))
    if(!all(metadataVars%in%coldatanames)) warning(" Can't find '",paste0(metadataVars[!metadataVars%in%coldatanames],collapse = ' & '),"' in coldata.\n '",paste0(metadataVars[metadataVars%in%coldatanames],collapse = ' & '), "' will be available for coloring plots in the app.")
    if(all(!metadataVars%in%coldatanames)) warning(" Can't find 'metadataVars' in coldata.\n Only 'partitionVars' will be available for coloring plots in the app.")
    coldatanames <- coldatanames[coldatanames%in%c("nCounts", "nFeatures", partitionVars, metadataVars)]
    colData(xx.sce) <- colData(xx.sce)[,coldatanames]
    # transform to character to factors to be able to plot in shiny app
    colsK <- sapply(colData(xx.sce), function(x){(is.character(x))})
    if(any(colsK)){
      colData(xx.sce)[,colsK] <- lapply(colData(xx.sce)[,colsK,drop=F], function(x){
      x <- as.factor(x)
      x <- droplevels(x)
      x})
	  }
  }

  # Attaching SCE to output ----
  csceo[["SCE"]] <- xx.sce
  
  
  if(verbose) cat('Computing differential expression markers:\n')
  
  # sce markers ----
  if(verbose) cat('Computing cluster markers...')
  # If test.type = wilcox, calculate FDR with that test and extract logFC values from t test.
  if (paramFindMarkers$test.type == "wilcox"){
    sce.markers <- list()
    for(i in ttoFactors){
      tout <- scran::findMarkers(xx.sce, 
                                  assay.type = "logcounts",
                                  group = colData(xx.sce)[,i],
                                  test.type="t",
                                  direction="any",
                                  pval.type="all",
                                  log.p=T,full.stats=T)
      wout <- scran::findMarkers(xx.sce, 
                                  assay.type = "logcounts",
                                  group = colData(xx.sce)[,i],
                                  test.type="wilcox",
                                  direction="any",
                                  pval.type="all",
                                  log.p=T,full.stats=T)
      l = length(wout)
      for (ii in 1:l){
        wout[[ii]][rownames(tout[[ii]]),"summary.stats"] = tout[[ii]][,"summary.stats"]
        it = grep("stats.", names(tout[[ii]]))
        for (jj in it){
          wout[[ii]][[jj]][rownames(tout[[ii]]),"AUC"] = tout[[ii]][[jj]][,"logFC"]
          names(wout[[ii]][[jj]]) = names(tout[[ii]][[jj]]) 
        }
      }
      sce.markers[[i]] = wout
    }
  }else{
    sce.markers <- list()
    for(i in ttoFactors){
      sce.markers[[i]] <- scran::findMarkers(xx.sce, 
                                      assay.type = "logcounts",
                                      group = colData(xx.sce)[,i],
                                      test.type=paramFindMarkers$test.type,
                                      direction="any",
                                      pval.type="all",
                                      log.p=T,full.stats=T)
    }
  }
  
  if(verbose) cat(' Finished\n')

  # Attaching sce.markers to output
  csceo[["sce.markers"]] <- sce.markers
  
  
  # ldf functions ----
  numCores <- max(1, parallel::detectCores() - 2, na.rm = TRUE)
  if(require("doParallel", quietly = T)) doParallel::registerDoParallel(numCores)

  ldf <- list()
  for(i in ttoFactors){
    ldf[[i]] <- ldf_func(xx.sce, i, paramFindMarkers)
  }
  if(verbose) cat('Finished\n')
  
  
  # Attaching ldf to output ----
  csceo[["ldf"]] <- ldf
  
  # Adding description text to output ----
  csceo[["text"]] <- descriptionText
  
  # Subsampling require? ----
  nmaxcell = 50000
  if(ncol(xx.sce)>nmaxcell){
      csceo$CELLS2KEEP <- subsampling_func(xx.sce, cellsToKeep = cells2keep, nmaxcell = nmaxcell)
  } else {
    csceo$CELLS2KEEP <- "all"
  }



  gc() # free memory
  return(csceo)
}

#' @import SingleCellExperiment

# Reduced Dimensions ----
#' @keywords internal
#' @noRd
applyReducedDim <- function(sce, reddimstocalculate, chosen.hvgs, nPCs, assayname, prefix.name="SCX_",verbose=TRUE){
  namepca <- paste0(prefix.name,"PCA")
  if("PCA"%in%reddimstocalculate){
    if(verbose) cat("\t PCA...")
    set.seed(12534)
    sce <- scater::runPCA(sce, subset_row=chosen.hvgs, ncomponents=nPCs, name=namepca, exprs_values=assayname)
    if(verbose) cat(' Finished','\n')
  }
  if("TSNE"%in%reddimstocalculate){
    if(verbose) cat("\t TSNE...")
    set.seed(1111011)
    sce <- scater::runTSNE(sce,dimred=namepca,n_dimred=20,ncomponents=3,name=paste0(prefix.name,"TSNE"),exprs_values=assayname)
    if(verbose) cat(' Finished','\n')
  }
  if("UMAP"%in%reddimstocalculate){
    if(verbose) cat("\t UMAP...")
    set.seed(1111011)
    sce <- scater::runUMAP(sce,dimred=namepca,n_dimred=20,ncomponents=3,name=paste0(prefix.name,"UMAP"),exprs_values=assayname)
    if(verbose) cat(' Finished','\n')
  }
  if("TSNE2D"%in%reddimstocalculate){
    if(verbose) cat("\t TSNE2D...")
    set.seed(1111011)
    sce <- scater::runTSNE(sce,dimred=namepca,name=paste0(prefix.name,"TSNE2D"),exprs_values=assayname)
    if(verbose) cat(' Finished','\n')
  }
  if("UMAP2D"%in%reddimstocalculate){
    if(verbose) cat("\t UMAP2D...")
    set.seed(1111011)
    sce <- scater::runUMAP(sce,dimred=namepca,name=paste0(prefix.name,"UMAP2D"),exprs_values=assayname)
    if(verbose) cat(' Finished','\n')
  }
  
  return(sce)
}

# ldf func for cluster markers ----
# returns: markers, robustness, correlation with a binary vector "turned on" in that cluster
#' @keywords internal
#' @noRd
ldf_func <- function(sce, partition, paramFindMarkers, minSize=50){

  cat(partition, ":\n", sep = "")
  
  # calculate lfmrk ----
  lfmrk <- list()
  lfmrk[[paramFindMarkers$pval.type]]    <- scran::findMarkers(sce,
                                                               assay.type = "logcounts",
                                                               groups=colData(sce)[,partition],
                                                               test.type=paramFindMarkers$test.type,
                                                               direction=paramFindMarkers$direction,
                                                               pval.type=paramFindMarkers$pval.type,
                                                               full.stats=TRUE
                                                               )
  
  
  # calculate ldf ----
  ldf_t <-list()
  
  lab   <- colData(sce)[,partition]
  tt   <- table(lab)
  labOK <- tt > minSize
  labOK <- names(labOK)[labOK]
  
  lfmrk <- lapply(lfmrk,function(x){x[labOK]})
  
  if(length(lfmrk[[1]])>0){ #At least one group to calculate
    for(ic in seq_along(lfmrk[[1]])){ #acá se podría hacer una paralelizacion
      coi <- names(lfmrk[[1]])[ic]
      cat('\t', coi,'- ')
      if(paramFindMarkers$pval.type=="any"){
        u <- rownames(lfmrk[[1]][[coi]])[lfmrk[[1]][[coi]][,'FDR']<0.05 & 
                                                lfmrk[[1]][[coi]][,'Top']<=10]
      } else {
        u <- rownames(lfmrk[[1]][[coi]])[lfmrk[[1]][[coi]][,'FDR']<0.05]
      }
      
      # (6.1.1) Boxcor ----
      cat("Computing correlation\n")
      
      if(length(u) > 0){
        Z   <- assay(sce, "counts")[u,,drop=FALSE]  
        
        pattern <- rep(0,ncol(Z))
        pattern[colData(sce)[,partition]%in%coi] <- 1
        boxcor  <-apply(Z,1,
                        function(x){
                          return(cor(x,pattern))
                        }
                       )

        df <- data.frame(summary.stats=lfmrk[[1]][[coi]][u,'summary.stats'],
                         log.FDR=lfmrk[[1]][[coi]][u,'FDR'],
                         boxcor=boxcor
                         )
        df   <- df[with(df, order(summary.stats, -log.FDR, boxcor, decreasing = T)),]
        df[] <- apply(df, 2, function(x){as.numeric(formatC(x, digits = 3))})
      }
      else {
        df <- NULL
      }
      ldf_t[[coi]]<-df
    }
  }
  return(ldf_t)
}

# subsampling cells for visual purposes
#' @keywords internal
#' @noRd
subsampling_func = function(sce, cellsToKeep=NULL, nmaxcell=50000){

  if(is.null(cellsToKeep)){
      ccells2keep <- sample(colnames(sce), nmaxcell, replace=FALSE)
      ccells2keep <- match(ccells2keep, colnames(sce))
  }else{
      ccells2keep <- cells2keep[cells2keep%in%colnames(sce)]
      if(length(ccells2keep)<nmaxcell){
        naddcells <- nmaxcell-length(cells2keep)
        ccells2keep <- c(ccells2keep, sample(colnames(sce)[!colnames(sce)%in%ccells2keep],naddcells,replace=FALSE))
        ccells2keep <- match(ccells2keep, colnames(sce))
      }
  }
  return(ccells2keep)
}

