# Managing input
#' @keywords internal
#' @noRd
inputManager <- function(xx, metadata){
  
  if(class(xx)[1]=="character"){
    xx <- DropletUtils::read10xCounts(samples = xx)
    colnames(xx) <- xx$Barcode
  }
  
  ## Check for colnames and rownames in the object
  if(is.null(rownames(xx)) | is.null(colnames(xx))){
    stop("Object must have rownames 'gene_id' and colnames 'barcodes'")
  } else {
  rnames <- rownames(xx)
    ##Check empty names
    if(length(which(rnames==""))>0){
      message("Some genes have empty names and will be excluded")
      rnames <- rnames[!rnames==""]
    }
    ##Check repeated names
    tt <- table(rnames)
    repnames <- names(tt[tt>1])
    if(length(repnames)>0){
      message("Some rownames are repeated and will be excluded: ", paste0(repnames, collapse = ', '))
      rnames <- rnames[!rnames%in%repnames]
    }
    xx <- xx[rnames,]
    
  ##Renaming if colnames are repeated
  cnames <- colnames(xx)
    tt <- table(cnames)
    repnames <- names(tt[tt>1])
    if(length(repnames)>0){
      message("Some colnames are repeated and will be renamed sequentially adding a '-' prior to the number")
      for(cname in repnames){
        ii <- which(cnames%in%cname)
        cnames[ii] <- paste0(cname,"-",seq_along(ii))
      }
      colnames(xx) <- cnames
    }
  }
  
  ## Before anything, stop if there is at least one cell not included in metadata
  if(!is.null(metadata)){
    if(!all(colnames(xx) %in% rownames(metadata))){
    	stop("Some cells in data are not present in metadata")
    }
  }
  return(xx)
}

## Creating SCE
#' @keywords internal
#' @noRd
sceConverter <- function(csceo, metadata, verbose){
  if(verbose) message('Creating SCE object... ', appendLF = F)
  ## Seurat to SCE object
  if(class(csceo$SCE)[1]=="Seurat"){
    ##Changing assay.name parameters because as.SingleCellExperiment fills 'counts' and 'logcounts' assays
    message("\nSeurat object is detected, 'assay.name.raw' & 'assay.name.normalization' will be set to default 'counts' & 'logcounts'")
    csceo$usage$assay.name.raw <- assay.name.raw <- "counts"
    csceo$usage$assay.name.normalization <- assay.name.normalization <- "logcounts"
    seurat.nms <- names(csceo$SCE@assays)
    ##Converts seurat object to sce object
    xx.sce <- Seurat::as.SingleCellExperiment(csceo$SCE)
  
    # if Seurat object contains only one matrix, as.SingleCellExperiment duplicates it as both 'counts' and 'logcounts'
    if(all(assay(xx.sce, "counts")==assay(xx.sce, "logcounts"))){
    	if(any(grep("logcount", seurat.nms))){
    		assay(xx.sce, "counts") <- NULL
    	}else{
    		assay(xx.sce, "logcounts") <- NULL
    	}
     }
    ##Converting seurat feature's metadata to rowData() of the sce object
    if(nrow(csceo$SCE@assays$RNA@meta.features)>0){
      rowData(xx.sce) <- csceo$SCE@assays$RNA@meta.features
    }
    ##If seurat object has variable features calculated and `chosen.hvg`=NULL then it will set to the variable features list
    if(length(Seurat::VariableFeatures(csceo$SCE))>0){
      if(is.null(csceo$call$chosen.hvg)){
        csceo$usage$chosen.hvg <- Seurat::VariableFeatures(csceo$SCE)
      }
    }
    
    for(nm in assayNames(xx.sce)){ # make every assay sparse
    	assay(xx.sce, nm) <- as(assay(xx.sce, nm),"dgCMatrix")
    }
  
  ## Matrix to SCE object
  } else if (class(csceo$SCE)[1] %in% c("dgCMatrix", "Matrix", "matrix")){
   xx.sce <- SingleCellExperiment(list(counts=as(csceo$SCE,"dgCMatrix")))
  
  ## xx as SCE object
  } else if (class(csceo$SCE)[1]=="SingleCellExperiment"){
      xx.sce <- csceo$SCE
      for(nm in assayNames(xx.sce)){ # make every assay sparse
      	assay(xx.sce, nm) <- as(assay(xx.sce, nm),"dgCMatrix")
      }
  } else { ##Stop if `xx` is not an object of the class Seurat, SCE or Matrix
  stop("xx must be an object of the class Seurat, SingleCellExperiment or Matrix")
  }
  
  ## Add metadata to sce if present
  if(!is.null(metadata)){
    ##Adding metadata to colData() of sce oject (colnames were checked before creating sce object)
    colData(xx.sce) <- cbind(colData(xx.sce), metadata[colnames(xx.sce),,drop=F])
  }
  
  if(verbose) message('Finished\nAnalyzing partitions... ', appendLF = F)
  return(xx.sce)
}

# partitionVars to factors
#' @keywords internal
#' @noRd
defPartitions <- function(csceo, verbose){
  ##Check for repeated partitions
  csceo$usage$partitionVars <- unique(csceo$call$partitionVars)
  
  ##Checking if partitions from 'partitionVars' are present in the metadata
  if(!all(csceo$usage$partitionVars %in% names(colData(csceo$SCE)))){ # Note: if partitionVars is NULL, condition is FALSE
    warning("At least one element from partitionVars is not present in metadata")
  }
  ## Check factors in colData
  tfs <- setNames(nm=csceo$usage$partitionVars,object = csceo$usage$partitionVars %in% names(colData(csceo$SCE)))
  
  if(any(tfs)){ # Note: if partitionVars is NULL, condition is FALSE
    ## Check if columns have more than one level
    vld <- sapply(colData(csceo$SCE)[,csceo$usage$partitionVars[tfs],drop=F],function(x){length(unique(x))>1})
    if(!any(vld)){ # if none of the partitionVars in colData have more than 1 level, quick cluster
  	  message("\npartitionVars in colData only have one level. A quick clusterization will be computed")
    }
    tfs[names(vld)] <-  tfs[names(vld)] & vld
    } else { # if partitionVars is NULL or not present in colData, quick cluster
        message("\npartitionVars unspecified or not found in colData. A quick clusterization will be computed")
    }
  
  ttoFactors <- names(tfs)[tfs] # partitions used for DEGs and markers
  ##If no partition passes the above filters
  if(length(ttoFactors) == 0){ # when partitionVars is NULL, not present in colData, or only have 1 level
    ##A quick clusterization will be computed in order to calculate gene markers and DEGs
    ttoFactors <- "scx.clust"
    colData(csceo$SCE)[ttoFactors] <- NA
  }else{
    ##Transforming partitions to factors
    colData(csceo$SCE)[ttoFactors] <- lapply(colData(csceo$SCE)[ttoFactors], as.factor)
  }
  
  # Checking number of levels of ttoFactors for calculations
  # When `calcAllPartitions=FALSE`:
  #   If none of the partitions in `partitionVars` have less than 31 levels, the function will ask user if continuing or using scx.clust
  #   If there is at least one partition with less than 31 levels, the function will continue by computing gene markers and DEGs only for those partitions while ignoring the rest.
  # When `calcAllPartitions=TRUE`:
  #   Compute gene markers and differentially expressed genes (DEGs) for all partitions in 'partitionVars' regardless of the maximum number of levels
  
  if(!csceo$usage$calcAllPartitions){
    allToFactors <- sapply(ttoFactors, function(x){length(unique(colData(csceo$SCE)[,x]))})>30
    if(all(allToFactors)){
  	  message("\n", paste0(names(allToFactors)[allToFactors], collapse =  ' & '), " have more than 30 levels.\nComputing markers and DEGs in this case could be very time-consuming, do you wish to proceed?")
  	  user <- readline("If not, a quick clusterization will be considered instead (y/n): ")
    	# if user input does not include the letter 'y', it is interpreted as 'no' and considers a quick cluster
    	if(length(grep("y", user, ignore.case = T))==0) ttoFactors <- "scx.clust"
      } else if(any(allToFactors)) {
          message("\n", paste0(paste0(names(allToFactors)[allToFactors], collapse =  ' & ')," have more than 30 levels. They wont be used to compute markers and DEGs. If you want to compute it anyway set 'calcAllPartitions' as TRUE"))
          ttoFactors <- ttoFactors[!allToFactors]
    }
  }
  if(verbose) message('Finished')
  
  csceo$usage$partitionVars <- ttoFactors
  if("scx.clust" %in% csceo$usage$partitionVars) csceo$usage$scx.clust = TRUE # if a quick clusterization is requierd or requested
  if("")
  
  return(csceo)
}

# QC
# Calculate the number of counts and features per cell
#' @keywords internal
#' @noRd
qcMetrics <- function(csceo, verbose){
  if(verbose) message('Computing QC metrics... ', appendLF = F)
  if(!csceo$usage$assay.name.raw %in% assayNames(csceo$SCE)){
    if(csceo$usage$assay.name.normalization %in% assayNames(csceo$SCE)){
      warning("Assay ", csceo$usage$assay.name.raw, " not found in SCE object")
      ## If there is no raw assay but the SCE object has a normalized assay then set `nCounts,nFeatures = NA`
      csceo$SCE$nCounts   <- NA
      csceo$SCE$nFeatures <- NA
    } else {
      ## If there is no raw or normalized assay available then the function will end
      stop("Assays ", paste(csceo$usage$assay.name.raw, csceo$usage$assay.name.normalization, sep = ' & '), " not found in SCE object")
    }
  } else {
      ## If there is a raw assay in the SCE object, calculate the number of counts
    csceo$SCE$nCounts <- colSums(assay(csceo$SCE, csceo$usage$assay.name.raw))
      ## 'Apply()' converts a sparse matrix to dense, forcing us to calculate the number of features per cell in an alternate manner
      if(class(assay(csceo$SCE, csceo$usage$assay.name.raw))[1]%in%c("dgCMatrix")){
        csceo$SCE$nFeatures <- diff(assay(csceo$SCE, csceo$usage$assay.name.raw)@p)
      } else {
        csceo$SCE$nFeatures <- colSums(assay(csceo$SCE, csceo$usage$assay.name.raw)>0)
      }
  }
  if(verbose) message('Finished')
  
  return(csceo)
}


# Normalization
# If there is no normalized assay in the SCE object:
#   Perform a 'Normalization by deconvolution' (proposed in the OSCA book [https://bioconductor.org/books/3.17/OSCA.basic/normalization.html#normalization-by-deconvolution])
#   First, calculate clusters using the Walktrap community detection algorithm for graph-based clustering with default parameters from `scran::quickCluster`.
#   The resulting clusters are stored in colData() as "scx.clust".
#   Then compute scale factors for the cells using the clusters.
#   Finally, calculate the lognormalized expression matrix by applying a log2 transformation to the product of the raw matrix and scale factors, with the addition of 1.
# If a normalized assay exists and "scx.clust" is included in the `partitionVars`, the function described before will be applied to compute the clusters, which are stored in colData()
#' @keywords internal
#' @noRd
expNormalization <- function(csceo, verbose){
  if(!csceo$usage$assay.name.normalization %in% assayNames(csceo$SCE)){
    if(csceo$usage$assay.name.raw %in% assayNames(csceo$SCE)){
      if(verbose) message('Computing normalization... ', appendLF = F)
      set.seed(123457)
      clust <- scran::quickCluster(csceo$SCE, assay.type = csceo$usage$assay.name.raw)
      csceo$SCE <- scran::computeSumFactors(csceo$SCE,cluster=clust,min.mean=0.1, assay.type = csceo$usage$assay.name.raw)
      csceo$SCE <- scater::logNormCounts(csceo$SCE, assay.type = csceo$usage$assay.name.raw, name="logcounts")
      if(verbose) message('Finished')
    }
  } 
  
  return(csceo)
}


# HVGs
# If `chosen.hvg` are not specified calculate the `nHVGs` most variable genes with biological component > 0
# This is computed with `scran::modelGeneVar` which calculates the variance and mean of the lognormalized expression values.
# By fitting a trend of the variance against the mean, a biological component of variation for each gene 
# can be assigned as the residual from the trend.
#' @keywords internal
#' @noRd
defHVG <- function(csceo, verbose){
    if(verbose) message('Computing HVGs... ', appendLF = F)
    mgv <- scran::modelGeneVar(csceo$SCE,span=.8, assay.type = csceo$usage$assay.name.normalization)
    rowData(csceo$SCE) <- cbind(rowData(csceo$SCE), hvg.mvBio=mgv$bio)
    chosen.hvg <- rank(-rowData(csceo$SCE)$hvg.mvBio) <= csceo$usage$nHVGs & rowData(csceo$SCE)$hvg.mvBio>0
    chosen.hvg <- rownames(csceo$SCE)[chosen.hvg]
    csceo$usage$chosen.hvg <- chosen.hvg
    
    if(verbose) message('Finished')
    
    return(csceo)
}


#' @keywords internal
#' @noRd
scXcluster <- function(csceo, blusparam, verbose){
  if(csceo$usage$scx.clust){
    if(verbose) message('Computing clusters... ', appendLF = F)
    if(is.null(blusparam)){ # if BLUSPARAM=NULL then it will use k and method from scx.graph.k/.method
      blusparam=bluster::NNGraphParam(k = csceo$usage$scx.graph.k, cluster.fun=csceo$usage$scx.graph.cluster.method)
    }
    if("SCX_PCA" %in% reducedDimNames(csceo$SCE)){
      csceo$SCE$scx.clust <- scran::clusterCells(csceo$SCE, use.dimred = "SCX_PCA", BLUSPARAM = blusparam)
    } else {
      csceo$SCE$scx.clust <- scran::clusterCells(csceo$SCE, assay.type = csceo$usage$assay.name.normalization, BLUSPARAM = blusparam)
    }
    if(verbose) message(' Finished\n')
  }
    
  return(csceo)
}

# Logcounts Normalized
#' @keywords internal
#' @noRd
lcNorm <- function(csceo, verbose){
  ##Changing name for normalized matrix to be input in shiny app
  if(csceo$usage$assay.name.normalization!="logcounts" & "logcounts" %in% assayNames(csceo$SCE)){
    assay(csceo$SCE, "logcounts") <- NULL
  }
  assayNames(csceo$SCE)[which(assayNames(csceo$SCE)==csceo$usage$assay.name.normalization)] <- "logcounts"
  ##Compute a row-normalization of the lognormalized expression matrix to be able to compare between gene expression profiles
  ##The row expression values are divided by their maximum value
  if(!"logcounts.norm" %in% assayNames(csceo$SCE)){
    if(verbose) message('Computing normalized logcounts... ', appendLF = F)
    sparse_mat <- as(assay(csceo$SCE, "logcounts"), "sparseMatrix")
    row_maxs <- sparseRowMax(sparse_mat)
    maxdiag <- Diagonal(x = 1/as.vector(row_maxs))
    scaled_sparse <- maxdiag %*% sparse_mat
    rownames(scaled_sparse) <- rownames(sparse_mat)
      
    ##Store the row-normalized expression values in `logcounts.norm` assay
    assays(csceo$SCE)$logcounts.norm <- scaled_sparse
    if(verbose) message('Finished')
  }
  
  return(csceo)
}

#' @keywords internal
#' @noRd
colPlotShiny <- function(csceo, verbose){
  # Transform columns to factors to be available for coloring plots in app
  if(!is.null(csceo$usage$metadataVars)){
  # Keep only names in `partitionVars` and `metadataVars` which are in colData
    coldatanames <- names(colData(csceo$SCE))
    if(!all(csceo$usage$metadataVars%in%coldatanames)) warning("Can't find ",paste0(csceo$usage$metadataVars[!csceo$usage$metadataVars%in%coldatanames],collapse = ' & ')," in data")
    
  # Subsetting SCE object
    colData(csceo$SCE) <- colData(csceo$SCE)[, coldatanames%in%c("nCounts", "nFeatures", "scx.clust", csceo$usage$partitionVars, csceo$usage$metadataVars)]
    
    colsK <- sapply(colData(csceo$SCE), function(x){(is.character(x) | is.numeric(x))})
    colsK[c("nCounts", "nFeatures")] <- F # added in QC and shouldn't be considered in any case
    if(any(colsK)){
      colData(csceo$SCE)[,colsK] <- lapply(colData(csceo$SCE)[,colsK,drop=F], function(x){
      x <- as.factor(x)
      x <- droplevels(x)
      x})
    }
  if(verbose) message(paste0(names(colData(csceo$SCE))[sapply(colData(csceo$SCE), function(x){(is.factor(x))})], collapse = ' & ')," will be available for coloring")
  } else {
  # Keep names in `partitionVars` which are in colData and any other column (character or numeric) with less than or equal to 30 levels
  colsK <- sapply(colData(csceo$SCE), function(x){(is.character(x) | is.numeric(x)) & length(unique(x)) <= 30})
  colsK[c("nCounts", "nFeatures")] <- F # added in QC and shouldn't be considered in any case
  if(any(colsK)){
  	colData(csceo$SCE)[,colsK] <- lapply(colData(csceo$SCE)[,colsK,drop=F], function(x){
  		x <- as.factor(x)
  		x <- droplevels(x)
  		x
  	})
  }
  }
  
  return(csceo)
}

# DEGs
#' @keywords internal
#' @noRd
degs <- function(csceo, verbose){
  if(verbose) message('Computing differential expression markers... ', appendLF = F)
  ##If test.type = wilcox, calculate FDR with that test and extract logFC values from t test.
  if (csceo$usage$paramFindMarkers$test.type == "wilcox"){
    sce.degs <- list()
    for(i in c(csceo$usage$partitionVars,"scx.clust")[c(csceo$usage$partitionVars,"scx.clust")%in%names(colData(csceo$SCE))]){
      tout <- scran::findMarkers(csceo$SCE, 
                                 assay.type = "logcounts",
                                 group = colData(csceo$SCE)[,i],
                                 test.type="t",
                                 direction="any",
                                 pval.type="all",
                                 log.p=T,full.stats=T,
                                 BPPARAM=csceo$usage$BPPARAM)
      wout <- scran::findMarkers(csceo$SCE, 
                                 assay.type = "logcounts",
                                 group = colData(csceo$SCE)[,i],
                                 test.type="wilcox",
                                 direction="any",
                                 pval.type="all",
                                 log.p=T,full.stats=T,
                                 BPPARAM=csceo$usage$BPPARAM)
      l = length(wout)
      for (ii in 1:l){
        wout[[ii]][rownames(tout[[ii]]),"summary.stats"] = tout[[ii]][,"summary.stats"]
        it = grep("stats.", names(tout[[ii]]))
        for (jj in it){
          wout[[ii]][[jj]][rownames(tout[[ii]]),"AUC"] = tout[[ii]][[jj]][,"logFC"]
          names(wout[[ii]][[jj]]) = names(tout[[ii]][[jj]]) 
        }
      }
      sce.degs[[i]] = wout
    }
  }else{
    sce.degs <- list()
    for(i in c(csceo$usage$partitionVars,"scx.clust")[c(csceo$usage$partitionVars,"scx.clust")%in%names(colData(csceo$SCE))]){
      sce.degs[[i]] <- scran::findMarkers(csceo$SCE,
                                          assay.type = "logcounts",
                                          group = colData(csceo$SCE)[,i],
                                          test.type=csceo$usage$paramFindMarkers$test.type,
                                          direction="any",
                                          pval.type="all",
                                          log.p=T,full.stats=T,
                                          BPPARAM = csceo$usage$BPPARAM)
    }
  }
  
  if(verbose) message('Finished\nComputing cluster markers...')
  return(sce.degs)
}

# Markers
# `sce.markers` is a list of lists of dataframes for every cluster in all partitions in `partitionVars` used by the shiny app to be able to identify gene markers for clusters.
# Every element in the list is a list of data The list is constructed using the output of `scran::findMarkers' with parameters specified by the user in `paramFindMarkers'.
# Only genes with FDR<0.05 are selected for each cluster, and the boxcor is calculated for those genes.
# Boxcor:
#   The boxcor is the correlation between a gene's expression vector (logcounts) and a binary vector, where only the cells from the selected cluster
#   mark 1 while the rest of the cells mark 0.
#' @keywords internal
#' @noRd
markers <- function(csceo, verbose){
  
  # if marker list
  if(!is.null(csceo$usage$markerList)){
    partitions <- unique(csceo$usage$markerList$Partition)
    if(!all(partitions %in% csceo$usage$partitionVars)){
      message(paste0(paste(partitions[!partitions %in% csceo$usage$partitionVars], collapse = ", ")," not present in 'partitionVars' and will be ignored.\n"), appendLF = F)
    }
  }
  
  # run marker_fun
  sce.markers <- list()
  for(i in c(csceo$usage$partitionVars,"scx.clust")[c(csceo$usage$partitionVars,"scx.clust")%in%names(colData(csceo$SCE))]){
    if(!is.null(csceo$usage$markerList)){
      if(i %in% partitions){
        sce.markers[[i]] <- markers_func(csceo$SCE, i, csceo$usage$markerList[csceo$usage$markerList$Partition==i,2:3], csceo$usage$paramFindMarkers, bpparam = csceo$usage$BPPARAM, minsize = csceo$usage$minSize, verbose = verbose)
      } else {
        sce.markers[[i]] <- markers_func(csceo$SCE, i, markerList = NULL, csceo$usage$paramFindMarkers, bpparam = csceo$usage$BPPARAM, minsize = csceo$usage$minSize, verbose = verbose)
      }
    } else { 
        sce.markers[[i]] <- markers_func(csceo$SCE, i, markerList = NULL, csceo$usage$paramFindMarkers, bpparam = csceo$usage$BPPARAM, minsize = csceo$usage$minSize, verbose = verbose)
    }
  }
  if(verbose) message('Finished')
  
  return(sce.markers)
}

#' @import SingleCellExperiment

# Reduced Dimensions
#' @keywords internal
#' @noRd
applyReducedDim <- function(sce, reddimstocalculate, chosen.hvgs, nPCs, assayname, prefix.name="SCX_",verbose=TRUE){
  if(verbose) message('Computing the following reduced dims: ',paste(reddimstocalculate, collapse = ' '))
  namepca <- paste0(prefix.name,"PCA")
  if("PCA"%in%reddimstocalculate){
    if(verbose) message("\tPCA... ", appendLF = F)
    set.seed(12534)
    sce <- scater::runPCA(sce, subset_row=chosen.hvgs, ncomponents=nPCs, name=namepca, exprs_values=assayname)
    if(verbose) message('Finished')
  }
  if("TSNE"%in%reddimstocalculate){
    if(verbose) message("\tTSNE... ", appendLF = F)
    set.seed(1111011)
    sce <- scater::runTSNE(sce,dimred=namepca,n_dimred=20,ncomponents=3,name=paste0(prefix.name,"TSNE"),exprs_values=assayname)
    if(verbose) message('Finished')
  }
  if("UMAP"%in%reddimstocalculate){
    if(verbose) message("\tUMAP... ", appendLF = F)
    set.seed(1111011)
    sce <- scater::runUMAP(sce,dimred=namepca,n_dimred=20,ncomponents=3,name=paste0(prefix.name,"UMAP"),exprs_values=assayname)
    if(verbose) message('Finished')
  }
  if("TSNE2D"%in%reddimstocalculate){
    if(verbose) message("\tTSNE2D... ", appendLF = F)
    set.seed(1111011)
    sce <- scater::runTSNE(sce,dimred=namepca,name=paste0(prefix.name,"TSNE2D"),exprs_values=assayname)
    if(verbose) message('Finished')
  }
  if("UMAP2D"%in%reddimstocalculate){
    if(verbose) message("\tUMAP2D... ", appendLF = F)
    set.seed(1111011)
    sce <- scater::runUMAP(sce,dimred=namepca,name=paste0(prefix.name,"UMAP2D"),exprs_values=assayname)
    if(verbose) message('Finished')
  }
  if(verbose) message('Finished')
  return(sce)
}

# markers func for cluster markers
# returns: markers, robustness, correlation with a binary vector "turned on" in that cluster
#' @keywords internal
#' @noRd
markers_func <- function(sce, partition, markerList, paramFindMarkers, bpparam, minsize=10, verbose = TRUE){ # previously ldf_func

  if(verbose) message(" ", partition, ":")
  
  # calculate lfmrk 
  lfmrk <- list()
  lfmrk[[paramFindMarkers$pval.type]]    <- scran::findMarkers(sce,
                                                               assay.type = "logcounts",
                                                               groups=colData(sce)[,partition],
                                                               test.type=paramFindMarkers$test.type,
                                                               direction=paramFindMarkers$direction,
                                                               pval.type=paramFindMarkers$pval.type,
                                                               full.stats=TRUE,
                                                               BPPARAM=bpparam
                                                               )
  
  
  # calculate sce.markers 
  scemarkers_t <-list()
  
  lab   <- colData(sce)[,partition]
  tt   <- table(lab)
  labOK <- tt > minsize
  labOK <- names(labOK)[labOK]
  
  lfmrk <- lapply(lfmrk,function(x){x[labOK]})
  
  if(length(lfmrk[[1]])>0){ #At least one group to calculate
    for(ic in seq_along(lfmrk[[1]])){ #acá se podría hacer una paralelizacion
      coi <- names(lfmrk[[1]])[ic]
      if(verbose) message('\t', coi,' - ', appendLF = F)
      
      if(paramFindMarkers$pval.type=="any"){
        u <- rownames(lfmrk[[1]][[coi]])[lfmrk[[1]][[coi]][,'FDR']<0.05 & 
                                                lfmrk[[1]][[coi]][,'Top']<=10]
      } else {
        u <- rownames(lfmrk[[1]][[coi]])[lfmrk[[1]][[coi]][,'FDR']<0.05]
      }
      if(!is.null(markerList)){
        u <- markerList[markerList$Cluster==coi, "Gene"]
        if(!all(u%in%rownames(sce))){
          message(paste0(paste(u[!u%in%rownames(sce)], collapse = ", "))," ignored because not in sce rownames.\n")
          u <- u[u%in%rownames(sce)]
        }
      }
      
      # (6.1.1) Boxcor 
      if(verbose) message("Computing correlation")
      
      if(length(u) > 0){
        Z   <- assay(sce, "logcounts")[u,,drop=FALSE]  
        
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
      scemarkers_t[[coi]]<-df
    }
  }
  return(scemarkers_t)
}

# subsampling cells for visual purposes
#' @keywords internal
#' @noRd
subsampling_func <- function(sce, cellsToKeep=NULL, nmaxcell=50000){

  if(is.null(cellsToKeep)){
      set.seed(123457) # seed
      ccells2keep <- sample(colnames(sce), nmaxcell, replace=FALSE)
      ccells2keep <- match(ccells2keep, colnames(sce))
  }else{
    ccells2keep <- cellsToKeep[cellsToKeep%in%colnames(sce)]
      if(length(ccells2keep)<nmaxcell){
        naddcells <- nmaxcell-length(ccells2keep)
        set.seed(123457)
        ccells2keep <- c(ccells2keep, sample(colnames(sce)[!colnames(sce)%in%ccells2keep],naddcells,replace=FALSE))
        ccells2keep <- match(ccells2keep, colnames(sce))
      }
  }
 
   return(ccells2keep)
}