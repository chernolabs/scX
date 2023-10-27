# Object to sce ----
#' SingleCellExperiment object ready for use in scX
#'
#' `createSCEobject()` returns a list which includes a SingleCellExperiment object with normalized expression, 
#'    gene markers and differntial expression analysis for each clusterization (at least one clusterization is computed), reduced dimensions for
#'		visualization, and any additional data provided. This list is used as input when launching
#'		the scX app.
#' 
#' @param xx Either a matrix with counts, or a SCE or Seurat object.
#' @param assay.name.raw Assay name for raw counts matrix if object is a SCE. Defaults to `counts`.
#' @param assay.name.normalization Assay name for normalized matrix if present in SCE. If not present,
#' 		it computes `logcounts` (default).
#' @param metadata Optional dataframe with cell metadata. Rownames in the metadata must include all cell names of `xx`.
#' @param partitionVars Optional metadata column names (or `colData`) to use for gene markers and differential
#'    expression analysis. If `NULL` (default),	a quick clusterization will be computed.
#' @param metadataVars Additional metadata column names (or `colData`) to use only for coloring in plots.
#'    If `NULL` (default), only `partitionVars` columns will be available for coloring plots.
#' @param chosen.hvg Optional list of Highly Variable Genes. NOTE: if `chosen.hvg=NULL` and `xx` is a Seurat object with computed VariableFeatures then 
#'    this parameter will be set to that list of genes.
#' @param nHVGs Number of Highly Variable Genes to use if `chosen.hvg=NULL`. Defaults to 3000.
#' @param nPCs Number of Principal Components to use in PCA. Defaults to 50.
#' @param calcRedDim Logical indicating whether to compute reduced dimensions (PCA, UMAP, TSNE, UMAP2D, 
#'		TSNE2D). Defaults to `TRUE`. If its set to `FALSE` but there is no 2D nor 3D or not at all reduced dimension calculated in the input:
#'    it calculates all. Instead if there is no 2D, it calculates PCA, UMAP2D, TSNE2D. And lastly if there is no 3D, it calculates PCA, UMAP, TSNE.
#' @param paramFindMarkers List of parameters to pass to `scran::findMarkers(...)` to compute marker 
#'		genes for clusters. Defaults to `list(test.type="wilcox", pval.type="all", direction="any")`.
#' @param calcAllPartitions Logical indicating whether to force the computation of markers and DEGs from
#'		the entire list of `partitionVars`. Defaults to `FALSE` and only partitions from `partitionVars` with less or equal to 30 levels will be account for
#'    markers and DEGs calculations.
#' @param cells2keep List of names for cells to keep when subsampling data. NOTE: Subsampling is only 
#'		used in large datasets for visual purposes, and it does not affect computations. Only 50k cells will be used for visualization in the app and
#'    their indexs are stored in the `CELLS2KEEP` element of the output list.
#' @param descriptionText Optional short description of the object being analized to be displayed in the `scX` app. This can help when 
#'		working with multiple tabs.
#' @param verbose Logical for step by step status while function is running. Defaults to `TRUE`.
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
  
  csceo <- list() # Creating list object to be returned by this function

  # Managing input ----
  # First we will ensure that there is a partition to be calculated

  ##Check if there are colnames and rownames in the objects 
  if(is.null(rownames(xx)) | is.null(colnames(xx))){
    stop('Missing row names or column names.')
  }
  ##Check for repeated partitions
  partitionVars <- unique(partitionVars)
  
  ##If there is not an specified 'partitionVars' then it sets to "scx.clust"
  if(is.null(partitionVars)){
    warning('No partitionVars specified, a quick clusterization will be computed.')
    partitionVars <- "scx.clust"
  }

  ##Also if it isn't a Seurat or SCE object, and metadata is null (metadata could be inside the object) the 'partitionVars' becomes "scx.clust"
  if(is.null(metadata) & class(xx)[1]!="Seurat" & class(xx)[1]!="SingleCellExperiment"){
    warning('No metadata specified, a quick clusterization will be computed.')
    partitionVars <- "scx.clust"
  }
  if(verbose) cat('Creating SCE object...')
  
  ## Seurat to SCE object ----
  ##Check if the input class is a Seurat object
  if(class(xx)[1]=="Seurat"){
    ##Changing assay.name parameters because as.SingleCellExperiment fills 'counts' and 'logcounts' assays
    warning("Seurat object is detected, 'assay.name.raw' & 'assay.name.normalization' will be set to default.")
    assay.name.raw <- "counts"
    assay.name.normalization <- "logcounts"
    ##Converts seurat object to sce object
    xx.sce <- Seurat::as.SingleCellExperiment(xx) 
    ##Checking if all partitions from 'partitionVars' are converted into the colData() of the sce object
    if((all(partitionVars!="scx.clust")) & (!all(partitionVars %in% names(colData(xx.sce))))){
      warning('at least one partition is not present in metadata')
    }
    ##Converting seurat feature's metadata to rowData() of the sce object
    if(nrow(xx@assays$RNA@meta.features)>0){
      rowData(xx.sce) <- xx@assays$RNA@meta.features
    }
    ##If seurat object has variable features calculated and `chosen.hvg`=NULL then it will set to the variable features list
    if(length(Seurat::VariableFeatures(xx))>0){
      if(is.null(chosen.hvg)){
        chosen.hvg <- Seurat::VariableFeatures(xx)
      }
    }
  } 
  ## Matrix+Metadata to SCE object ----
  ##Check if the input class is a dense or sparse matrix object
  else if (class(xx)[1] %in% c("dgCMatrix", "Matrix", "matrix")){
    ##The matrix must contain rownames and colnames for the function to be able to identify genes and cells
    if(is.null(rownames(xx)) | is.null(colnames(xx))){
      stop('Matrix must have rownames "gene_id" and colnames "barcodes"')
    }
    ##Creating a SCE object from the matrix
    xx.sce <- SingleCellExperiment(list(counts=xx))
    ##Checking if metadata was specified
    if(!is.null(metadata)){
      if(all(colnames(xx.sce) %in% rownames(metadata))){
        ##Assing colData() of the sce oject to be the metadata if all the cells have metadata
        colData(xx.sce) <- cbind(colData(xx.sce), metadata[colnames(xx.sce),,drop=F])
        ##Checking if partitions from 'partitionVars' are present in the metadata
        if((!("scx.clust" %in% partitionVars)) & (!all(partitionVars %in% names(colData(xx.sce))))){
          warning('at least one partition is not present in metadata')
        }
      ##Stop if there is at least one cell it is not included in the metadata
      } else {
        stop('Some cells in metadata are not present in the matrix colnames.\n')
        }
    }
  } 
  ## xx as SCE object ----
  ##Check if the input class is a SCE object
  else if (class(xx)[1]=="SingleCellExperiment"){
    xx.sce <- xx
    ##Checking if metadata was specified
    if(!is.null(metadata)){
      if(all(colnames(xx.sce) %in% rownames(metadata))){
        ##Adding to colData() of the sce oject the metadata if all the cells have metadata
        colData(xx.sce) <- cbind(colData(xx.sce), metadata[colnames(xx.sce),,drop=F])
        ##Checking if partitions from 'partitionVars' are present in the metadata
        if((partitionVars!="scx.clust") & (!all(partitionVars %in% names(colData(xx.sce))))){
          warning('at least one partition is not present in metadata')
        }
        ##Stop if there is at least one cell it is not included in the metadata
      } else {
        stop('Some cells in metadata are not present in the SCE colnames.\n')
      }
    } 
  ##Stop if the `xx` is not an object of the class Seurat, SCE or Matrix
  } else {
    stop('xx must be an object of the class Seurat, SingleCellExperiment or Matrix')
  }
  if(verbose) cat(' Finished\n')
  

  # partitionVars to factors ----
  # Setting which partitions will be transformed to factors

  if(verbose) cat('Changing factors from partitionVars...')
  ##Check the factors that are in the colData.
  tfs <- setNames(nm=partitionVars,object = partitionVars %in% names(colData(xx.sce)))
  ##Check if columns have more than one level.
  if(sum(tfs) > 0){
    vld <- sapply(colData(xx.sce)[,partitionVars[tfs],drop=F],function(x){length(unique(x))>1})   
    tfs[names(vld)] <-  tfs[names(vld)] & vld
  }
  
  ttoFactors <- names(tfs)[tfs]
  ##Warning message if a partition is not found in colData or if it has only one level.
  if(sum(!tfs)>0) warning("Can't find ",paste0(names(tfs)[!tfs],collapse = ' & ')," in metadata or they have only one level")
  ##If no partition passes the above filters
  if(length(ttoFactors) == 0){
    warning("No partition passed the controls, a quick clusterization will be computed.")
    ##A quick clusterization will be computed in order to calculate gene markers and DEGs
    ttoFactors <- "scx.clust"
  }else{
    ##Transforming partitions to factors (in case a partition is numeric or character)
	  colData(xx.sce)[ttoFactors] <- lapply(colData(xx.sce)[ttoFactors], as.factor)
  }
  

  # metadata ----
  #Transform all the character columns to factor to be able to be selected in the shinyApp if they have less or equal to 30 levels. 

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
  # When `calcAllPartitions=FALSE`:
  #   If none of the partitions in `partitionVars` have less than 31 levels, the function will stop.
  #   However, if there is at least one partition with less than 31 levels, the function will continue by computing gene markers and DEGs only for those partitions while ignoring the rest.
  # When `calcAllPartitions=TRUE`:
  #   Compute gene markers and differentially expressed genes (DEGs) for all partitions in 'partitionVars' regardless of the maximum number of levels

  if(!calcAllPartitions){
    allToFactors <- sapply(ttoFactors, function(x){length(unique(colData(xx.sce)[,x]))})>30
    if(all(allToFactors)){
      stop(paste0(paste0(names(allToFactors)[allToFactors], collapse =  ' & ')," has more than 30 levels. If you want to compute it anyway set 'calcAllPartitions' as TRUE"))
    } else if(any(allToFactors)) {
      warning(paste0(paste0(names(allToFactors)[allToFactors], collapse =  ' & ')," has more than 30 levels. They wont be used to compute markers and DEGs. If you want to compute it anyway set 'calcAllPartitions' as TRUE"))
      ttoFactors <- ttoFactors[!allToFactors]
    }
  }
  

  # QC ----
  # Caluclate the number of counts and features per cell

  if(verbose) cat('Computing QC metrics...')
  if(!assay.name.raw %in% names(assays(xx.sce))){
    if(assay.name.normalization %in% names(assays(xx.sce))){
      warning(paste0('Assays ',assay.name.raw,' not found in SCE object'))
      ## If there is no raw assay but the SCE object has a normalized assay then set `nCounts,nFeatures = NA`
      xx.sce$nCounts   <- NA
      xx.sce$nFeatures <- NA
    } else {
      ## If there is no raw or normalized assay available then the function will end
      stop(paste0('Assay ',paste(assay.name.raw, assay.name.normalization, sep = ' & '),' not found in SCE object'))
    }
  } else {
      ## If there is a raw assay in the SCE object, calculate the number of counts
      xx.sce$nCounts <- colSums(assay(xx.sce, assay.name.raw))
      ## 'Apply()' converts a sparse matrix to dense, forcing us to calculate the number of features per cell in an alternate manner
      if(class(assay(xx.sce, assay.name.raw))[1]%in%c("dgCMatrix")){
          xx.sce$nFeatures <- diff(assay(xx.sce, assay.name.raw)@p)
      } else {
          xx.sce$nFeatures <- apply(assay(xx.sce, assay.name.raw),2,function(x){sum(x>0)})
      }
  }
  if(verbose) cat(' Finished\n')
  
  
  # Normalization ----
  # If there is no normalized assay in the SCE object:
  #   Perform a 'Normalization by deconvolution' (proposed in the OSCA book [https://bioconductor.org/books/3.17/OSCA.basic/normalization.html#normalization-by-deconvolution])
  #   First, calculate clusters using the Walktrap community detection algorithm for graph-based clustering with default parameters from `scran::quickCluster`.
  #   The resulting clusters are stored in colData() as "scx.clust".
  #   Then compute scale factors for the cells using the clusters.
  #   Finally, calculate the lognormalized expression matrix by applying a log2 transformation to the product of the raw matrix and scale factors, with the addition of 1.
  # If a normalized assay exists and "scx.clust" is included in the `partitionVars`, the function described before will be applied to compute the clusters, which are stored in colData()

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
  } else if ( "scx.clust" %in% partitionVars ) {
      if(verbose) cat('Computing clusters...')
      clust <- scran::quickCluster(xx.sce, assay.type = assay.name.normalization)
      xx.sce$scx.clust <- clust
      if(verbose) cat(' Finished\n')
  }
  
  
  # HVGs ----
  # If `chosen.hvg` are not specified calculate the `nHVGs` most variable genes with biological component > 0
  # This is computed with `scran::modelGeneVar` which calculates the variance and mean of the lognormalized expression values.
  # By fitting a trend of the variance against the mean, a biological component of variation for each gene 
  # can be assigned as the residual from the trend.

  if(is.null(chosen.hvg)){
    if(verbose) cat('Computing HVGs...')
    mgv <- scran::modelGeneVar(xx.sce,span=.8, assay.type = assay.name.normalization)
    rowData(xx.sce) <- cbind(rowData(xx.sce), hvg.mvBio=mgv$bio)
    chosen.hvg <- rank(-rowData(xx.sce)$hvg.mvBio) <= nHVGs & rowData(xx.sce)$hvg.mvBio>0
    chosen.hvg <- rownames(xx.sce)[chosen.hvg]
    if(verbose) cat(' Finished\n')
  }
  
  
  # Reduced dimensions ----
  # If user specifies no calcRedDim but:
  #   there is no dimRed calculated -> it calculates all
  #   there is no 2D or 3D redDim caclulated -> it calculates all
  #   there is no 2D -> it calculates only 2D
  #   there is no 3D -> it calculates only 3D
  # PCA is calculated with `scater::runPCA()` using the `chosen.hvg` and retaining the first `nPCs` components
  # TSNE and UMAP are calculated with `scater::runTSNE()` and `scater::runUMAP()` functions using the PCA matrix

  runDim <- c("PCA", "TSNE", "UMAP", "TSNE2D", "UMAP2D")
  if(!calcRedDim){
    if(length(reducedDimNames(xx.sce))<1 | all(!(sapply(reducedDims(xx.sce),ncol) %in% c(2,3)))){
      calcRedDim <- TRUE
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
  

  # Logcounts Normalized ----
  
  ##Changing name for normalized matrix to be input in shiny app
  if(assay.name.normalization!="logcounts" & "logcounts" %in% assayNames(xx.sce)){
    assay(xx.sce, "logcounts") <- NULL
  }
  assayNames(xx.sce)[which(assayNames(xx.sce)==assay.name.normalization)] <- "logcounts"
  ##Compute a row-normalization of the lognormalized expression matrix to be able to compare between gene expression profiles
  ##The row expression values are divided by their maximum value
  if(verbose) cat('Computing logcounts normalized...')
  if(!"logcounts.norm" %in% names(assays(xx.sce))){
		sparse_mat <- as(assay(xx.sce, "logcounts"), "sparseMatrix")
		row_maxs <- qlcMatrix::rowMax(sparse_mat)
		maxdiag <- Diagonal(x = 1/as.vector(row_maxs))
		scaled_sparse <- maxdiag %*% sparse_mat
		rownames(scaled_sparse) <- rownames(sparse_mat)
    ##Store the row-normalized expression values in `logcounts.norm` assay
		assays(xx.sce)$logcounts.norm <- scaled_sparse
  }
  if(verbose) cat(' Finished\n')
  

  # Subsetting SCE object ----
  # Keep only names in `partitionVars` and `metadataVars` which are in colData to be used for colouring plots
  # Transform character to factors to be able to plot in shiny app

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
  
  
  # DEGs ----
  # List of dataframes used by the shiny app to be able to plot volcano graphs of DEGs between clusters.
  # The list contains the output of `scran::findMarkers' by setting the direction to "any" and pval.type to "all", 
  # the test.type is specified by the user in `paramFindMarkers', with the exception that if the user 
  # had specified "wilcox" as test.type, the logFC values would be extracted from the function with test.type="t".
  
  if(verbose) cat('Computing differential expression markers:\n')
  if(verbose) cat('Computing cluster markers...')
  ##If test.type = wilcox, calculate FDR with that test and extract logFC values from t test.
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


  # Attaching sce.markers to output ----
  csceo[["sce.markers"]] <- sce.markers
  
  
  # Markers ----
  # `ldf` is a list of lists of dataframes for every cluster in all partitions in `partitionVars` used by the shiny app to be able to identify gene markers for clusters.
  # Every element in the list is a list of data The list is constructed using the output of `scran::findMarkers' with parameters specified by the user in `paramFindMarkers'.
  # Only genes with FDR<0.05 are selected for each cluster, and the boxcor is calculated for those genes.
  # Boxcor:
  #   The boxcor is the correlation between a gene's expression vector (logcounts) and a binary vector, where only the cells from the selected cluster
  #   mark 1 while the rest of the cells mark 0.

  numCores <- max(1, parallel::detectCores() - 2, na.rm = TRUE)
  if(require("doParallel", quietly = T)) doParallel::registerDoParallel(numCores)

  ldf <- list()
  for(i in ttoFactors){
    ldf[[i]] <- ldf_func(xx.sce, i, paramFindMarkers)
  }
  if(verbose) cat('Finished\n')
  
  
  # Attaching ldf to output ----

  csceo[["ldf"]] <- ldf
  

  # Attaching text to output ----
  # `descriptionText` is useful when working with multiple tabs because it is displayed in the browser hosting the app

  csceo[["text"]] <- descriptionText
  

  # Attaching CELLS2KEEP to output ----
  # If the SCE object contains over 50,000 cells, a random sample of 50,000 cells will be chosen for visualization in the application. 
  # Please note that all calculations are already completed, and this step is solely for promoting smooth and efficient visualization.

  nmaxcell = 50000
  if(ncol(xx.sce)>nmaxcell){
      if(verbose) cat('Subsampling sce object, it exceeds 50k cells\n')
      csceo$CELLS2KEEP <- subsampling_func(xx.sce, cellsToKeep = cells2keep, nmaxcell = nmaxcell)
      if(verbose) cat('Finished\n')
  } else {
    csceo$CELLS2KEEP <- "all"
  }


  ##Free memory
  gc() 
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