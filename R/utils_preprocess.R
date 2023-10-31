#' Creates an object ready for use in scX
#'
#' @description 
#' \code{\link{createSCEobject}} creates the input object (as a \linkS4class{List}) for the function \code{\link{launch_scX}}. This list includes a \linkS4class{SingleCellExperiment} object with normalized expression, 
#' reduced dimensions for visualization, and any additional data provided. Also, it includes gene markers and differential expression analysis for each specified partition (if no partition is selected, a quick clusterization is computed)
#' 
#' @param xx Either a numeric matrix-like object of the raw count matrix (with genes as rows and cells as columns), a \linkS4class{SingleCellExperiment} object or \linkS4class{Seurat} object.
#' @param assay.name.raw Assay name for raw counts matrix if object is a \linkS4class{SingleCellExperiment}. Defaults to \code{counts}.
#' @param assay.name.normalization Assay name for normalized matrix if present in \linkS4class{SingleCellExperiment}. If not present,
#' 		it computes \code{logcounts} (default).
#' @param metadata Optional \linkS4class{DataFrame} with cell metadata. Rownames in the metadata must include all cell names of \code{xx}.
#' @param partitionVars Optional metadata column names (or \code{colData}) to use for gene markers and differential
#'    expression analysis. If \code{NULL} (default),	a quick clusterization using \code{\link[scran]{quickCluster}} from \pkg{scran} package will be computed.
#' @param metadataVars Additional metadata column names (or \code{colData})) to use only for coloring in plots.
#'    If \code{NULL} (default), only \code{partitionVars} columns will be available for coloring plots.
#' @param chosen.hvg Optional list of Highly Variable Genes. NOTE: if \code{chosen.hvg=NULL} and \code{xx} is a \linkS4class{Seurat} object with computed \code{\link[Seurat]{VariableFeatures}} then 
#'    this parameter will be set to that list of genes.
#' @param nHVGs Number of Highly Variable Genes to use if \code{chosen.hvg=NULL}. Defaults to 3000.
#' @param nPCs Number of Principal Components to use in PCA. Defaults to 50.
#' @param calcRedDim Logical indicating whether to compute reduced dimensions (PCA, UMAP, TSNE, UMAP2D, 
#'		TSNE2D). Defaults to \code{TRUE}. If its set to \code{FALSE} but there is no 2D nor 3D or not at all reduced dimension calculated in the input:
#'    it calculates all. Instead if there is no 2D, it calculates PCA, UMAP2D, TSNE2D. And lastly if there is no 3D, it calculates PCA, UMAP, TSNE.
#' @param paramFindMarkers List of parameters to pass to \code{\link[scran]{findMarkers}} to compute marker 
#'		genes for clusters. Defaults to \code{list(test.type="wilcox", pval.type="all", direction="any")}.
#' @param calcAllPartitions Logical indicating whether to force the computation of markers and DEGs from
#'		the entire list of \code{partitionVars}. Defaults to \code{FALSE} and only partitions from \code{partitionVars} with less or equal to 30 levels will be account for
#'    markers and DEGs calculations.
#' @param cells2keep List of names for cells to keep when subsampling data. NOTE: Subsampling is only 
#'		used in large datasets for visual purposes, and it does not affect computations. Only 50k cells will be used for visualization in the app and
#'    their indexs are stored in the \code{CELLS2KEEP} element of the output list.
#' @param descriptionText Optional short description of the object being analized to be displayed in the \code{scX} app. This can help when 
#'		working with multiple tabs.
#' @param verbose Logical for step by step status while function is running. Defaults to \code{TRUE}.

#' @details 
#' This function handle the basic preprocessing steps for sc/sn-RNAseq experiments data taking advantage of the utility power of \pkg{scran}, \pkg{scater} and \pkg{SingleCellExperiment} packages.
#' The steps are: Converting the input to \linkS4class{SingleCellExperiment} object, calculate some QC metrics, normalization of the expression matrix, find the most highly variable genes,
#' compute embeddings for various dimensionality reduction techniques, compute if not specified a clustering of the cells for marker genes and DEGs analysis,
#' subsampling the \linkS4class{SingleCellExperiment} object to reduce waiting times between plots in the \pkg{scX} app. Below you can find more details for these steps that are recommended in 
#' \href{https://bioconductor.org/books/release/OSCA/.}{OSCA} book.

#' @section Partitions:
#' It is extremely recommended that a curated partition of the data be present in the \linkS4class{SingleCellExperiment} object or in the \code{metadata}.
#' Marker genes and differential expression analysis will be compute for partitions that are passed through \code{partitionVars} and have less than 31 groups.
#' If user wants run the calculations even if some partitions have more than 30 levels must set \code{calcAllPartitions} to \code{TRUE}. 
#' It may be worth coloring some partitions in plots but without having to compute marker genes and DEGs analyses for them. Those partitions
#' had to be passed through \code{metadataVars}.
#' If \code{partitionVars} and \code{metadataVars} are both \code{NULL} a quick clustering (see "Normalization") is used for the marker genes and DEGs analysises.

#' @section Normalization:
#' If there is no normalized assay in the \linkS4class{SingleCellExperiment} object or \linkS4class{Seurat} object:
#' Perform a \href{https://bioconductor.org/books/3.17/OSCA.basic/normalization.html#normalization-by-deconvolution}{Normalization by deconvolution} proposed in the OSCA book.
#' First, calculate clusters using the Walktrap community detection algorithm for graph-based clustering with default parameters from \code{\link[scran]{quickCluster}}.
#' The resulting clusters are stored in \code{colData} as \code{"scx.clust"}.
#' Then compute scale factors for the cells using those clusters.
#' Finally, calculate the lognormalized expression matrix by applying a log2 transformation to the product of the raw matrix and scale factors with pseudocounts of 1.
#' If a normalized assay exists and "scx.clust" is included in the \code{partitionVars}, the function described before will be applied to compute the clusters, which are stored in \code{colData}.

#' @section Highly Variable Genes:
#' If \code{chosen.hvg} is not specified, it will calculate the top \code{nHVGs} most variable genes with biological component > 0
#' This is computed with \code{\link[scran]{modelGeneVar}} which calculates the variance and mean of the lognormalized expression values. 
#' By fitting a trend of the variance against the mean, a biological component of variation for each gene can be assigned as the residual from the trend.
#' See appropiate documentation for more details.

#' @section Reduced Dimensional Analysis Techniques:
#' \pkg{scX} is a tool that helps visualizing the data properly, so it runs some dimensionality reduction analysis technique (PCA, UMAP2D, TSNE2D, UMAP3D, TSNE3D)
#' If \code{xx} already has dimensionality reduction embeddings calculated, \code{calcRedDim} can be set to \code{FALSE} but:
#' - if there is no dimRed calculated in the \linkS4class{SingleCellExperiment} object -> it calculates all
#' - if there is no 2D or 3D redDim caclulated -> it calculates all
#' - if there is no 2D -> it calculates only 2D
#' - if there is no 3D -> it calculates only 3D
#' Principal Component Analysis is calculated with \code{\link[scater]{runPCA}} using the \code{chosen.hvg} and retaining the first \code{nPCs} components,
#' for the normalized expression matrix.
#' TSNE and UMAP are calculated with \code{\link[scater]{runTSNE}} and \code{\link[scater]{runUMAP}} functions using the PCA matrix.

#' @section Marker Genes Analysis:
#' For all partitions in \code{partitionVars} \link{createSCEobject} compute statistics for identify marker genes for every cluster
#' using \code{\link[scran]{findMarkers}} with parameters specified by the user in \code{paramFindMarkers}.
#' Only genes with FDR<0.05 are selected for each cluster, and the boxcor is calculated for those genes.
#' \describe{
#' \item{\code{boxcor}:}{The boxcor is the correlation between a gene's expression vector (logcounts) and a binary vector, where only the cells from the selected cluster
#' mark 1 while the rest of the cells mark 0.}
#' }

#' @section Differential Expression Analysis:
#' For all partitions in \code{partitionVars} \link{createSCEobject} compute statistics for identify DEGs between clusters using \code{\link[scran]{findMarkers}} by setting 
#' \code{direction="any"} and \code{pval.type="all"}, the test.type is specified by the user in \code{paramFindMarkers}, with the exception that if the user 
#' had specified "wilcox" as \code{test.type}, the logFC values would be extracted from the function with \code{test.type="t"}.

#' @section Subsampling cells:
#' If the \linkS4class{SingleCellExperiment} object contains over 50k cells, a random sample of 50k cells will be chosen for visualization in the application.
#' Cell names that the user wants to keep in the visualizations can be passed to \code{cells2keep}.
#' Please note that all calculations are already completed, and this step is solely for promoting smooth and efficient visualization.

#' @return 
#' A named \linkS4class{List} that serves as input for \code{\link{launch_scX}} which contains the fields:
#' \describe{
#' \item{\code{SCE}:}{\linkS4class{SingleCellExperiment} object with a computed normalized expression,
#' dimensional reduction embedding (PCA, UMAP, TSNE, in 2D & 3D) calculated using the list of \code{chosen.hvg} if not \code{NULL} or 
#' the top \code{nHVGs}. \code{colData} contains the \code{partitionVArs} and \code{metadataVars} if they where specified, if not only
#' a quick clusterization will be available for preliminar analyisis of the data.}
#' \item{\code{sce.degs}:}{A named \linkS4class{List} of \linkS4class{DataFrame}s where each DataFrame contains the consolidated marker statistics
#' for each gene (row) for the cluster of the same name. Statistics are computed using \code{\link[scran]{findMarkers}} and user can choose \code{test.type}
#' parameters to pass to that function. See \code{\link[scran]{combineMarkers}} for the details of how these dataframes are confectioner.}
#' \item{\code{sce.markers}:}{A \linkS4class{List} of named \linkS4class{List}s of \linkS4class{DataFrame}s. Each one corresponds to the marker genes of every cluster in a partition 
#' (names of the nested lists). \code{summary.stats:} AUC if \code{test.type=="wilcox"} and -log.FC for \code{test.type=="t"} or \code{test.type=="binom"}.
#' \code{log.FDR:} -Log.FDR of the most appropriate inter-cluster comparison according to the selected p-value type. See \code{\link[scran]{findMarkers}} for the details of how these metrics are computed.
#' \code{boxcor:} Correlation scores between the normalized gene expression profiles and a binary vector of cells, in which cells of the selected cluster have a value of 1.}
#' \item{(optional) \code{text}:}{String if \code{descriptionText} not \code{NULL}.}
#' \item{\code{CELLS2KEEP}:}{Numeric, the indices of the selected cells chosen for visualization in the \code{scX} application. 
#' (Alternatively) Character, if 'all' no subsampling will be performed for visualization purposes.}
#' }

#' @author 
#' Tomás Vega Waichman, Maria Luz Vercesi, Ariel A. Berardino, Maximiliano S. Beckel, Chernomoretz Lab and Collaborators

#' @seealso 
#' Related functions from \pkg{scran} and \pkg{scater} packages as suggested in \href{https://bioconductor.org/books/release/OSCA/.}{OSCA} book for:
#' \itemize{
#' \item{Preprocessing steps: }{\code{\link[scran]{quickCluster}}, \code{\link[scran]{computeSumFactors}} and \code{\link[scater]{logNormCounts}}.}
#' \item{HVGs: }{\code{\link[scran]{modelGeneVar}}.}
#' \item{Dimensionality Reduction Techniques: }{\code{\link[scater]{runPCA}}, \code{\link[scater]{runTSNE}} and \code{\link[scater]{runUMAP}}.}
#' \item{Marker genes and DEGs analyses: }{\code{\link[scran]{findMarkers}}, \code{\link[scran]{combineMarkers}}.}
#' }

#' @examples 
#' # Quick start guide example:
#' library(scX)
#' cseo <- createSCEobject(xx = sce, 
#'                         partitionVars = "inferred_cell_type", 
#'                         metadataVars = c("source_name", "age", "sex", "strain", "treatment", "pseudotime"),
#'                         descriptionText = "Quick Start Guide")
#' launch_scX(cseo)

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
  ## Matrix+Metadata to SCE object ----
  ##Check if the input class is a dense or sparse matrix object
  } else if (class(xx)[1] %in% c("dgCMatrix", "Matrix", "matrix")){
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
    ## xx as SCE object ----
    ##Check if the input class is a SCE object
  } else if (class(xx)[1]=="SingleCellExperiment"){
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
    colData(xx.sce)[ttoFactors] <- NA
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
    sce.degs <- list()
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
      sce.degs[[i]] = wout
    }
  }else{
    sce.degs <- list()
    for(i in ttoFactors){
      sce.degs[[i]] <- scran::findMarkers(xx.sce, 
                                      assay.type = "logcounts",
                                      group = colData(xx.sce)[,i],
                                      test.type=paramFindMarkers$test.type,
                                      direction="any",
                                      pval.type="all",
                                      log.p=T,full.stats=T)
    }
  }
  
  if(verbose) cat(' Finished\n')


  # Attaching sce.degs to output ----
  csceo[["sce.degs"]] <- sce.degs
  
  
  # Markers ----
  # `sce.markers` is a list of lists of dataframes for every cluster in all partitions in `partitionVars` used by the shiny app to be able to identify gene markers for clusters.
  # Every element in the list is a list of data The list is constructed using the output of `scran::findMarkers' with parameters specified by the user in `paramFindMarkers'.
  # Only genes with FDR<0.05 are selected for each cluster, and the boxcor is calculated for those genes.
  # Boxcor:
  #   The boxcor is the correlation between a gene's expression vector (logcounts) and a binary vector, where only the cells from the selected cluster
  #   mark 1 while the rest of the cells mark 0.

  numCores <- max(1, parallel::detectCores() - 2, na.rm = TRUE)
  if(require("doParallel", quietly = T)) doParallel::registerDoParallel(numCores)

  sce.markers <- list()
  for(i in ttoFactors){
    sce.markers[[i]] <- markers_func(xx.sce, i, paramFindMarkers)
  }
  if(verbose) cat('Finished\n')
  
  
  # Attaching sce.markers to output ----

  csceo[["sce.markers"]] <- sce.markers
  

  # Attaching text to output ----
  # `descriptionText` is useful when working with multiple tabs because it is displayed in the browser hosting the app

  if(!is.null(descriptionText)){
    if(class(descriptionText)=="character"){
      csceo[["text"]] <- descriptionText
    } else {
      warning("'descriptionText' is not a character vector.. setting 'text' to NULL")
    }
  }
  

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

# markers func for cluster markers ----
# returns: markers, robustness, correlation with a binary vector "turned on" in that cluster
#' @keywords internal
#' @noRd
markers_func <- function(sce, partition, paramFindMarkers, minSize=50){ # previously ldf_func

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
  
  
  # calculate sce.markers ----
  scemarkers_t <-list()
  
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
      scemarkers_t[[coi]]<-df
    }
  }
  return(scemarkers_t)
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
