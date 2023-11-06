#' Creates an object ready for use in the scX app
#' 
#' @description 
#' \code{\link{createSCEobject}} creates the input object (a \linkS4class{List}) for the function \code{\link{launch_scX}}. This list includes a \linkS4class{SingleCellExperiment} object with normalized expression, 
#' reduced dimensions for visualization, and any additional data provided. Also, this object includes identified gene markers and differential expression analysis for each user-specified partition (if no partition was selected, a quick clusterization is computed and considered)
#' 
#' @param xx Either a numeric count matrix object (with genes as rows and cells as columns), a \linkS4class{SingleCellExperiment} object or a \linkS4class{Seurat} object.
#' @param assay.name.raw Assay name for the raw counts matrix if the input object is a \linkS4class{SingleCellExperiment}. Defaults to \code{counts}.
#' @param assay.name.normalization Assay name for the normalized matrix if present in the \linkS4class{SingleCellExperiment}. If not present, it computes \code{logcounts} by default.
#' @param metadata (Optional) A \linkS4class{DataFrame} containing cell metadata. The row names of the metadata data frame must include all cell names that appear as columns in \code{xx}.
#' @param partitionVars (Optional) Names of metadata or \code{colData} columns to be used for gene markers and differential
#' expression analysis. If set to \code{NULL} (default), a quick clustering step will be performed using \code{\link[scran]{quickCluster}} from the \pkg{scran} package.
#' @param metadataVars (Optional) Names of additional metadata or \code{colData} columns to be used for coloring in plots.
#' If set to \code{NULL} (default), only \code{partitionVars} columns will be available for coloring plots.
#' @param chosen.hvg (Optional) A list of Highly Variable Genes. NOTE: If \code{chosen.hvg=NULL} and \code{xx} is a \linkS4class{Seurat} object with computed \code{\link[Seurat]{VariableFeatures}}, then this parameter will be set to that list of genes.
#' @param nHVGs Number of Highly Variable Genes to use if \code{chosen.hvg=NULL}. Defaults to 3000.
#' @param nPCs Number of Principal Components to use in PCA. Defaults to 50.
#' @param calcRedDim Logical. Indicates whether to compute reduced dimensions (PCA, UMAP, t-SNE, UMAP2D, t-SNE2D). Defaults to \code{TRUE}. (Note: If set to \code{FALSE}, but there are no 2D and >3D reduced dimensions provided in the input, PCA embedding will be estimated anyway).
#' @param paramFindMarkers A named list of parameters to pass to \code{\link[scran]{findMarkers}} to compute cluster gene markers. 
#' Defaults to \code{list(test.type="wilcox", pval.type="all", direction="up")}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating whether and how parallelization should be performed across genes in the \code{\link[scran]{findMarkers}} function.
#' @param minSize Numeric. The minimum cluster size for calculating gene marker statistics.
#' @param calcAllPartitions Logical. Defaults to \code{FALSE}, which means that only partitions from \code{partitionVars} with 30 or fewer levels will be considered for marker and DEG calculations. If set to TRUE, it forces the computation of markers and DEGs for the entire list of \code{partitionVars}.
#' @param cells2keep (Optional) A list of cell names to keep in case of subsampling. NOTE: Subsampling is only activated for visualization purposes in the case of large datasets; it is not used for computations. 50k cells will be used for visualization in the app, and their indexes are stored in the \code{CELLS2KEEP} element of the CSEO object.
#' @param nSubCells Numeric. The maximum number of cells to select for subsampling the data set.
#' @param descriptionText (Optional) A short description of the object being analyzed. This text will be displayed in the Summary module of the \code{scX} app.
#' @param verbose Logical. Indicates whether to show step-by-step status while the function is running. Defaults to \code{TRUE}.

#' @details 
#' This function handles the basic preprocessing steps for sc/sn-RNAseq data. It leverages functionality implemented in the \pkg{scran}, \pkg{scater}, and \pkg{SingleCellExperiment} packages.
#' The steps include: converting the input to a \linkS4class{SingleCellExperiment} object, estimating QC metrics, normalizing the expression matrix, identifying highly variable genes,
#' and computing embeddings for various dimensionality reduction techniques. If no cell partition is provided, a default clustering step is conducted to investigate marker genes and differential expression patterns.
#' For large datasets (with more than 50k cells), subsampling of the \linkS4class{SingleCellExperiment} object is considered to reduce waiting times in the \pkg{scX} app. A detailed discussion of the employed preprocessing pipeline can be found in the
#' \href{https://bioconductor.org/books/release/OSCA/.}{OSCA} book.

#' @section Partitions:
#' It is recommended to include a curated partition of the data in the input object or in the \code{metadata}.
#' Marker genes and differential expression analysis will be automatically computed for partitions with fewer than 31 levels specified as \code{partitionVars}.
#' If the user wishes to run the calculations for partitions with more than 30 levels, \code{calcAllPartitions} must be set to \code{TRUE} (Please note that this step could be very time-consuming).
#' You may want to color some partitions in plots without having to compute marker genes and DEGs analyses for them. To do this, pass those partitions through \code{metadataVars}.

#' @section Metadata:
#' Metadata is passed to \link{createSCEobject} as a \linkS4class{DataFrame} where the rows must include all the cells present in the \code{XX} input. Metadata columns represent cell covariates. All character or numeric covariates passed in \code{metadataVars} with fewer than 31 unique values are set to factors and are referred to as "Categories" in the scX app. Numeric covariates with more than 30 unique values are called "Fields."

#' @section Normalization:
#' If there is no normalized assay in the \linkS4class{SingleCellExperiment} object or \linkS4class{Seurat} input object
#' a \href{https://bioconductor.org/books/3.17/OSCA.basic/normalization.html#normalization-by-deconvolution}{Normalization by deconvolution} is performed as proposed in the OSCA book.
#' First, we calculate clusters using the Walktrap community detection algorithm for graph-based clustering (default parameters from \code{\link[scran]{quickCluster}}).
#' The resulting clusters are stored in \code{colData} as \code{"scx.clust"}.
#' Then we compute scale factors for the cells using those clusters.
#' Finally, we calculate the lognormalized expression matrix by applying a log2 transformation to the product of the raw matrix and scale factors considering a pseudocounts of 1.

#' @section Clustering:
#' If \code{partitionVars} is user-specified, those categories are used for clustering the data.
#' If \code{partitionVars} is \code{NULL}, the quick clustering step (see "Normalization") is used to identify markers and DEGs.
#' Additionally, if a normalized assay exists in the SCE object and "scx.clust" is included in the \code{partitionVars}, the previously described function will be applied to compute clustering, which is stored in \code{colData} as "scx.clust."

#' @section Highly Variable Genes:
#' If \code{chosen.hvg} is not specified, we will use \code{\link[scran]{modelGeneVar}} to calculate the variance and mean of the lognormalized expression values. 
#' By fitting a trend of the variance against the mean, a biological component of variation for each gene can be assigned as the residual from the trend (see \pkg{scran} documentation for more details).
#' We consider the top \code{nHVGs} most variable genes with a biological component greater than 0.

#' @section Reduced Dimensional Analysis Techniques:
#' \pkg{scX} is a tool that helps visualize the data properly by running several dimensionality reduction analysis techniques (PCA, UMAP2D, TSNE2D, UMAP3D, TSNE3D). To use the app, \code{xx} must have at least a 2D dimensional reduction embedding.
#' If \code{xx} already has dimensionality reduction embeddings calculated, you can set \code{calcRedDim} to \code{FALSE}. NOTE: If set to \code{FALSE}, but there are no dimension reductions calculated in the \linkS4class{SingleCellExperiment} object, or if they all have less than 4 dimensions with none of them having 2D, a PCA will be calculated.
#' Principal Component Analysis is calculated with \code{\link[scater]{runPCA}} using the \code{chosen.hvg} and retaining the first \code{nPCs} components for the normalized expression matrix.
#' TSNE and UMAP are calculated with \code{\link[scater]{runTSNE}} and \code{\link[scater]{runUMAP}} functions using the PCA matrix.

#' @section Marker Genes Analysis:
#' \link{createSCEobject} computes statistics to identify marker genes for every cluster
#' for all partitions in \code{partitionVars}  
#' using \code{\link[scran]{findMarkers}} functions (with parameters specified by the user in \code{paramFindMarkers}).
#' Only genes with FDR < 0.05 are selected for each cluster. The boxcor score is calculated for those genes as follows:
#' \describe{
#' \item{\code{boxcor}:}{The boxcor is the correlation between a gene's expression vector (logcounts) and a binary vector, where only the cells from the selected cluster
#' mark 1 while the rest of the cells mark 0.}
#' }

#' @section Differential Expression Analysis:
#' We use the \code{\link[scran]{findMarkers}} function to identify DEGs between clusters specified in \code{partitionVars}
#' (\code{direction="any"},\code{pval.type="all"}).
#' The test.type can be specified by the user in \code{paramFindMarkers$test.type}.

#' @section Subsampling cells:
#' If the \linkS4class{SingleCellExperiment} object contains over than 50k cells (this number can be changed through \code{nSubCells}), a random sample of 50k cells will be chosen for visualization purposes in the application.
#' Cell names that the user wants to keep in the visualizations can be specified in the \code{cells2keep} parameter.
#' Please note that it is solely for producing efficient visualizations.

#' @return 
#' A named \linkS4class{List} that serves as input for \code{\link{launch_scX}}, which contains the following fields:
#' \describe{
#' \item{\code{SCE}:}{A \linkS4class{SingleCellExperiment} object with computed normalized expression and dimensional reduction embeddings (PCA, UMAP, t-SNE, in 2D & 3D). These are calculated using the list of \code{chosen.hvg} if not \code{NULL}, or the top \code{nHVGs}. \code{colData} contains the \code{partitionVArs} and \code{metadataVars} if they were specified; otherwise, only a quick clusterization will be available for preliminary analysis of the data.}
#' \item{\code{sce.degs}:}{A named \linkS4class{List} of \linkS4class{DataFrame}s where each DataFrame contains the consolidated marker statistics for each gene (row) for the cluster of the same name. Statistics are computed using \code{\link[scran]{findMarkers}}, and the user can choose the \code{test.type} parameters to pass to that function. See \code{\link[scran]{combineMarkers}} for the details of how these dataframes are created.}
#' \item{\code{sce.markers}:}{A \linkS4class{List} of named \linkS4class{List}s of \linkS4class{DataFrame}s. Each one corresponds to the marker genes of every cluster in a partition (names of the nested lists). \code{summary.stats:} AUC if \code{test.type=="wilcox"} and -log.FC for \code{test.type=="t"} or \code{test.type=="binom"}. \code{log.FDR:} -Log.FDR of the most appropriate inter-cluster comparison according to the selected p-value type. See \code{\link[scran]{findMarkers}} for the details of how these metrics are computed. \code{boxcor:} Correlation scores between the normalized gene expression profiles and a binary vector of cells, in which cells of the selected cluster have a value of 1.}
#' \item{(optional) \code{text}:}{String if \code{descriptionText} not \code{NULL}.}
#' \item{\code{CELLS2KEEP}:}{Numeric, the indices of the selected cells chosen for visualization in the \code{scX} application. 
#' (Alternatively) Character, if 'all' no subsampling will be performed for visualization purposes.}
#' }

#' @author 
#' Tomás Vega Waichman, Maria Luz Vercesi, Ariel A. Berardino, Maximiliano S. Beckel, Chernomoretz Lab and Collaborators

#' @seealso 
#' Related functions from \pkg{scran} and \pkg{scater} packages as suggested in the \href{https://bioconductor.org/books/release/OSCA/.}{OSCA} book for:
#' \itemize{
#' \item{Preprocessing steps: }{\code{\link[scran]{quickCluster}}, \code{\link[scran]{computeSumFactors}} and \code{\link[scater]{logNormCounts}}.}
#' \item{HVGs: }{\code{\link[scran]{modelGeneVar}}.}
#' \item{Dimensionality Reduction Techniques: }{\code{\link[scater]{runPCA}}, \code{\link[scater]{runTSNE}} and \code{\link[scater]{runUMAP}}.}
#' \item{Marker genes and DEGs analyses: }{\code{\link[scran]{findMarkers}}, \code{\link[scran]{combineMarkers}}.}
#' }

#' @examples 
#' \dontrun{
#' # Quick start guide example:
#' library(scX)
#' cseo <- createSCEobject(xx = sce, 
#'                         partitionVars = "inferred_cell_type", 
#'                         metadataVars = c("source_name", "age", "sex", "strain", "treatment", "pseudotime"),
#'                         descriptionText = "Quick Start Guide")
#' launch_scX(cseo)
#' }

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
                            paramFindMarkers=list(test.type="wilcox", pval.type="all", direction="up"),
                            BPPARAM=BiocParallel::SerialParam(),
                            minSize=30,
                            calcAllPartitions=FALSE,
                            cells2keep=NULL,
                            nSubCells=50000,
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
    ##Transforming partitions to factors
	  colData(xx.sce)[ttoFactors] <- lapply(colData(xx.sce)[ttoFactors], as.factor)
  }
  

  # metadata ----
  #Transform all the character columns to factor to be able to be selected in the shinyApp if they have less or equal to 30 levels or they are characters. 

  cols <- sapply(colData(xx.sce), function(x){(is.character(x) | is.numeric(x)) & length(unique(x)) <= 30})
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
  # If calcRedDim = TRUE: it calculates PCA, TSNE, UMAP, TSNE2D, UMAP2D
  # If calcRedDim = FALSE: 
  #   if there is no dimRed calculated or none dimRed has more than 3 columns -> it calculates PCA
  #   except there is a dimRed 2D calculated
  # PCA is calculated with `scater::runPCA()` using the `chosen.hvg` and retaining the first `nPCs` components
  # TSNE and UMAP are calculated with `scater::runTSNE()` and `scater::runUMAP()` functions using the PCA matrix

  runDim <- c("PCA", "TSNE", "UMAP", "TSNE2D", "UMAP2D")
  if(!calcRedDim){
    if(length(reducedDimNames(xx.sce))<1){
      runDim <- c("PCA")
      calcRedDim = TRUE
    } else if( (all(sapply(reducedDims(xx.sce),ncol) < 4)) & (all(sapply(reducedDims(xx.sce),ncol) != 2))){
      runDim <- c("PCA")
      calcRedDim = TRUE
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
  
  if(!'test.type'%in%names(paramFindMarkers)) paramFindMarkers$test.type <- 'wilcox'
  if(!'pval.type'%in%names(paramFindMarkers)) paramFindMarkers$pval.type <- 'all'
  if(!'direction'%in%names(paramFindMarkers)) paramFindMarkers$direction <- 'up'
  
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
                                 log.p=T,full.stats=T,
                                 BPPARAM=BPPARAM)
      wout <- scran::findMarkers(xx.sce, 
                                 assay.type = "logcounts",
                                 group = colData(xx.sce)[,i],
                                 test.type="wilcox",
                                 direction="any",
                                 pval.type="all",
                                 log.p=T,full.stats=T,
                                 BPPARAM=BPPARAM)
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
                                          log.p=T,full.stats=T,
                                          BPPARAM = BPPARAM)
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
  sce.markers <- list()
  for(i in ttoFactors){
    sce.markers[[i]] <- markers_func(xx.sce, i, paramFindMarkers, bpparam = BPPARAM, minsize = minSize)
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

  if(ncol(xx.sce)>nSubCells){
      if(verbose) cat('Subsampling sce object, it exceeds 50k cells\n')
      csceo$CELLS2KEEP <- subsampling_func(xx.sce, cellsToKeep = cells2keep, nmaxcell = nSubCells)
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
markers_func <- function(sce, partition, paramFindMarkers, bpparam, minsize=10){ # previously ldf_func

  cat(partition, ":\n", sep = "")
  
  # calculate lfmrk ----
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
  
  
  # calculate sce.markers ----
  scemarkers_t <-list()
  
  lab   <- colData(sce)[,partition]
  tt   <- table(lab)
  labOK <- tt > minsize
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
