#' Creates an object ready for use in the scX app
#'
#' @description 
#' \code{\link{createSCEobject}} creates the input object (a \linkS4class{List}) for the function \code{\link{launch_scX}}. This list includes a \linkS4class{SingleCellExperiment} object with normalized expression, 
#' reduced dimensions for visualization, and any additional data provided. Also, this object includes identified gene markers and differential expression analysis for each user-specified partition (if no partition was selected, a quick clusterization is computed and considered)
#' 
#' @param xx Either a numeric count matrix object (with genes as rows and cells as columns), a \linkS4class{SingleCellExperiment} object, a \linkS4class{Seurat} object or the directory containing the matrix.mtx, features.tsv, barcodes.tsv provided by CellRanger.
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
#' @param markerList (Optional) A \linkS4class{DataFrame} with three columns named: 'Partition', 'Cluster', 'Gene'. It will be used to subset the cluster marker list shown in cluster markers in Markers tab in the scX app.
#' @param paramFindMarkers A named list of parameters to pass to \code{\link[scran]{findMarkers}} to compute cluster gene markers. 
#' Defaults to \code{list(test.type="wilcox", pval.type="all", direction="up")}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating whether and how parallelization should be performed across genes in the \code{\link[scran]{findMarkers}} function.
#' @param scx.clust Logical. Indicates whether to perform a community detection algorithm on the data. Defaults to \code{FALSE}. (Note: If set to \code{FALSE}, but there is no \code{partitionVars} to compute marker analysis, it will set to \code{TRUE}).
#' @param scx.graph.k An integer scalar specifying the number of nearest neighbors to consider during graph construction. Defaults to 20
#' @param scx.graph.cluster.method Function specifying the method to use to detect communities in the shared NN graph. Alternatively, this may be a string containing the suffix of any \pkg{igraph} community detection algorithm (for example: "walktrap" will use \code{\link[igraph]{cluster_walktrap}}). Defaults to "louvain".
#' @param BLUSPARAM A \linkS4class{BlusterParam} object specifying the clustering algorithm to use, defaults to a graph-based method. If not \code{NULL} this parameter will be used instead \code{scx.graph.k} and \code{scx.graph.cluster.method}.
#' @param minSize Numeric. The minimum cluster size for calculating gene marker statistics.
#' @param calcAllPartitions Logical. Defaults to \code{FALSE}, which means that only partitions from \code{partitionVars} with 30 or fewer levels will be considered for marker and DEG calculations. If set to TRUE, it forces the computation of markers and DEGs for the entire list of \code{partitionVars}.
#' @param cells2keep (Optional) A list of cell names to keep in case of subsampling. NOTE: Subsampling is only activated for visualization purposes in the case of large datasets; it is not used for computations. Only \code{nSubCells} cells will be used for visualization in the app, and their indexes are stored in the \code{CELLS2KEEP} element of the CSEO object.
#' @param nSubCells Numeric. The maximum number of cells for randomly subsampling the data set (just for visualization purposes).
#' @param descriptionText (Optional) A short description of the object being analyzed. This text will be displayed in the Summary module of the \code{scX} app.
#' @param verbose Logical. Indicates whether to show step-by-step status while the function is running. Defaults to \code{TRUE}.

#' @details 
#' This function handles the basic preprocessing steps for sc/sn-RNAseq data. It leverages functionality implemented in the \pkg{scran}, \pkg{scater}, and \pkg{SingleCellExperiment} packages.
#' The steps include: converting the input to a \linkS4class{SingleCellExperiment} object, estimating QC metrics, normalizing the expression matrix, identifying highly variable genes,
#' and computing embeddings for various dimensionality reduction techniques. If no cell partition is provided, a default clustering step is conducted to investigate marker genes and differential expression patterns.
#' For datasets with more than 50k cells, subsampling of the \linkS4class{SingleCellExperiment} object is suggested to improve smooth interactive visualizations in the \pkg{scX} app. A detailed discussion of the employed preprocessing pipeline can be found in the
#' \href{https://bioconductor.org/books/release/OSCA/.}{OSCA} book.

#' @section Partitions:
#' It is recommended to include a curated partition of the data in the input object or in the \code{metadata}.
#' Marker genes and differential expression analysis will be automatically computed for partitions with fewer than 31 levels specified as \code{partitionVars}.
#' If the user wishes to run the calculations for partitions with more than 30 levels, \code{calcAllPartitions} must be set to \code{TRUE} (Please note that this step could be very time-consuming).
#' You may want to color some partitions in plots without having to compute marker genes and DEGs analyses for them. To do this, pass those partitions as \code{metadataVars}.

#' @section Metadata:
#' Metadata is passed to \link{createSCEobject} as a \linkS4class{DataFrame} where the rows must include all the cells present in the \code{XX} input. Metadata columns represent cell covariates. All character or numeric covariates passed in \code{metadataVars} with fewer than 31 unique values are set to factors and are referred to as "Categories" in the scX app. Numeric covariates with more than 30 unique values are called "Fields."

#' @section Normalization:
#' If there is no normalized assay in the \linkS4class{SingleCellExperiment} object or \linkS4class{Seurat} input object
#' a \href{https://bioconductor.org/books/3.17/OSCA.basic/normalization.html#normalization-by-deconvolution}{Normalization by deconvolution} is performed as proposed in the OSCA book.
#' First, we calculate clusters using the Walktrap community detection algorithm for graph-based clustering (default parameters from \code{\link[scran]{quickCluster}}).
#' Then we compute scale factors for the cells using those clusters.
#' Finally, we calculate the lognormalized expression matrix by applying a log2 transformation to the product of the raw matrix and scale factors considering a pseudocounts of 1.

#' @section Highly Variable Genes:
#' If \code{chosen.hvg} is not specified, we will use \code{\link[scran]{modelGeneVar}} to calculate the variance and mean of the lognormalized expression values. 
#' By fitting a trend of the variance against the mean, a biological component of variation for each gene can be assigned as the residual from the trend (see \pkg{scran} documentation for more details).
#' We consider the top \code{nHVGs} most variable genes with a biological component greater than 0.

#' @section Reduced Dimensional Analysis Techniques:
#' \pkg{scX} is a tool that helps visualize the data properly by running several dimensionality reduction analysis techniques (PCA, UMAP2D, TSNE2D, UMAP3D, TSNE3D). To use the app, \code{xx} must have at least a 2D dimensional reduction embedding.
#' If \code{xx} already has dimensionality reduction embeddings calculated, you can set \code{calcRedDim} to \code{FALSE}. NOTE: If set to \code{FALSE}, but there are no dimension reductions calculated in the \linkS4class{SingleCellExperiment} object, or if they all have less than 4 dimensions with none of them having 2D, a PCA will be calculated.
#' Principal Component Analysis is calculated with \code{\link[scater]{runPCA}} using the \code{chosen.hvg} and retaining the first \code{nPCs} components for the normalized expression matrix.
#' TSNE and UMAP are calculated with \code{\link[scater]{runTSNE}} and \code{\link[scater]{runUMAP}} functions using the PCA matrix.

#' @section Clustering:
#' If \code{partitionVars} is user-specified, those categories are used for clustering the data.
#' If \code{partitionVars} is \code{NULL}, none of the user-specified categories pass the filter (see "Partitions" for details) or \code{scx.clust} is set to \code{TRUE} a clustering performed with a graph-based community detection algorithm will be included in partitionVars to identify markers and DEGs. By default it computes \code{\link[igraph]{cluster_louvain}} on a shared nearest neighbors graph with \code{k=20}. Both the clustering algorithm and the number of first neighbors can be changed using the \code{scx.graph.cluster.method} and \code{scx.graph.k} parameters. Also, the user can specify any \linkS4class{BlusterParam} object through \code{BBLUSPARAM} that will be passed to \code{\link[scran]{clusterCells}} to compute clusters.

#' @section Marker Genes Analysis:
#' \link{createSCEobject} computes statistics to identify marker genes for every cluster
#' for all partitions in \code{partitionVars}  
#' using \code{\link[scran]{findMarkers}} functions (with parameters specified by the user in \code{paramFindMarkers}).
#' Only genes with FDR < 0.05 are selected for each cluster. The boxcor score is calculated for those genes as follows:
#' \describe{
#' \item{\code{boxcor}:}{The boxcor is the correlation between a gene's expression vector (logcounts) and a binary vector, where only the cells from the selected cluster
#' mark 1 while the rest of the cells mark 0.}
#' }
#' The user can provide a \linkS4class{DataFrame} containing marker genes for each cluster in any partition through \code{markerList}. For each partition in \code{markerList} specified in \code{partitionVars}, the calculated marker statistics for each cluster will be subset to the genes in \code{markerList}. Note that the statistics were calculated using all genes.

#' @section Differential Expression Analysis:
#' We use the \code{\link[scran]{findMarkers}} function to identify DEGs between clusters specified in \code{partitionVars}
#' (\code{direction="any"},\code{pval.type="all"}).
#' The test.type can be specified by the user in \code{paramFindMarkers$test.type}.

#' @section Subsampling cells:
#' If the \linkS4class{SingleCellExperiment} object contains over \code{nSubCells} cells (50k by default), a random sample of that size will be chosen for visualization purposes in the application.
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
#' \item{\code{call}:}{A named list containing all parameters and their values passed to the \code{createSCEobject} call.}
#' \item{\code{usage}:}{A named list containing all parameters and their corresponding values that were finally utilized by \code{createSCEobject} during preprocessing.}

#' }

#' @author 
#' Tomás Vega Waichman, Maria Luz Vercesi, Ariel A. Berardino, Maximiliano S. Beckel, Chernomoretz Lab and Collaborators

#' @seealso 
#' Related functions from \pkg{scran} and \pkg{scater} packages as suggested in the \href{https://bioconductor.org/books/release/OSCA/.}{OSCA} book for:
#' \itemize{
#' \item{Preprocessing steps: }{\code{\link[scran]{quickCluster}}, \code{\link[scran]{computeSumFactors}} and \code{\link[scater]{logNormCounts}}.}
#' \item{HVGs: }{\code{\link[scran]{modelGeneVar}}.}
#' \item{Dimensionality Reduction Techniques: }{\code{\link[scater]{runPCA}}, \code{\link[scater]{runTSNE}} and \code{\link[scater]{runUMAP}}.}
#' \item{Clustering: }{\code{\link[scran]{clusterCells}}, \linkS4class{BlusterParam}}
#' \item{Marker genes and DEGs analyses: }{\code{\link[scran]{findMarkers}}, \code{\link[scran]{combineMarkers}}.}
#' }

#' @examples 
#' \dontrun{
#' # Quick start guide example:
#' library(scX)
#' cseo <- createSCEobject(xx = scXample, 
#'                         partitionVars = "inferred_cell_type", 
#'                         metadataVars = c("source_name", "age", "sex", "strain", "treatment", "pseudotime"),
#'                         descriptionText = "Quick Start Guide")
#' launch_scX(cseo)
#' }
#' 
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
                            markerList=NULL,
                            paramFindMarkers=list(test.type="wilcox", pval.type="all", direction="up"),
                            BPPARAM=BiocParallel::SerialParam(),
                            scx.clust=FALSE,
                            scx.graph.k=20,
                            scx.graph.cluster.method="louvain",
                            BLUSPARAM=NULL,
                            minSize=30,
                            calcAllPartitions=FALSE,
                            cells2keep=NULL,
                            nSubCells=50000,
                            descriptionText=NULL,
                            verbose=TRUE){
  
  csceo <- sapply(c("SCE", "sce.degs", "sce.markers", "CELLS2KEEP", "call", "usage"), function(x) NULL) # Creating list object to be returned by this function
  
  csceo$call <- list(
    assay.name.raw = assay.name.raw,
    assay.name.normalization = assay.name.normalization,
    partitionVars = partitionVars,
    metadataVars = metadataVars,
    chosen.hvg = chosen.hvg,
    nHVGs = nHVGs,
    nPCs = nPCs,
    calcRedDim = calcRedDim,
    markerList = markerList, 
    paramFindMarkers = paramFindMarkers,
    BPPARAM = BPPARAM,
    scx.clust=scx.clust,
    scx.graph.k=scx.graph.k,
    scx.graph.cluster.method=scx.graph.cluster.method,
    BLUSPARAM = BLUSPARAM,
    minSize = minSize,
    calcAllPartitions = calcAllPartitions,
    cells2keep = cells2keep,
    nSubCells = nSubCells,
    descriptionText = descriptionText)
  
  csceo$usage <- csceo$call
  
  on.exit(sink())
  logdir <- file.path(getwd(),"scXlogs") # subdir for every preprocess log
  dir.create(logdir, showWarnings = F)
  
  xxname <- deparse(substitute(xx))
  logfile <- paste0(logdir,.Platform$file.sep, xxname, "_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".log")
  sink(logfile)
  cat("scX preprocessing steps for object:", xxname, "\n*---------------------------------*\n")
  
  # Managing input ----
  csceo$SCE <- inputManager(xx, metadata)
  
  # input to SCE ----
  csceo$SCE <- sceConverter(csceo, metadata, verbose)
  
  ## partitionVars to factors ----
  csceo <- defPartitions(csceo, verbose)
  
  # QC ----
  csceo <- qcMetrics(csceo,verbose)
  
  # Normalization ----
  csceo <- expNormalization(csceo, verbose)
 
  # HVGs ----
  if(is.null(csceo$call$chosen.hvg)){
    csceo <- defHVG(csceo, verbose)
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
    if(length(reducedDimNames(csceo$SCE))<1){
      runDim <- c("PCA")
      calcRedDim = TRUE
    } else if( (all(sapply(reducedDims(csceo$SCE),ncol) < 4)) & (all(sapply(reducedDims(csceo$SCE),ncol) != 2))){
      runDim <- c("PCA")
      calcRedDim = TRUE
    }
  }
  
  if(calcRedDim){
    csceo$SCE <- applyReducedDim(csceo$SCE, 
                                 runDim, 
                                 csceo$usage$chosen.hvg, 
                                 csceo$usage$nPCs, 
                                 csceo$usage$assay.name.normalization, 
                                 prefix.name="SCX_", 
                                 verbose)
  }
  
  # clustering ----
  csceo <- scXcluster(csceo, csceo$usage$BLUSPARAM, verbose)
  
  # Logcounts Normalized ----
  csceo <- lcNorm(csceo, verbose)
  
  # Coloring plots ----
  csceo <- colPlotShiny(csceo, verbose)
  
  # DEGs ----
  if(!'test.type'%in%names(csceo$usage$paramFindMarkers)) csceo$usage$paramFindMarkers$test.type <- 'wilcox'
  if(!'pval.type'%in%names(csceo$usage$paramFindMarkers)) csceo$usage$paramFindMarkers$pval.type <- 'all'
  if(!'direction'%in%names(csceo$usage$paramFindMarkers)) csceo$usage$paramFindMarkers$direction <- 'up'
  
  csceo$sce.degs <- degs(csceo, verbose)
  
  # Markers ----
  csceo$sce.markers <- markers(csceo, verbose)

  
  # *Attaching text to output ----
  # `descriptionText` is useful when working with multiple tabs because it is displayed in the browser hosting the app
  
  if(!is.null(csceo$usage$descriptionText)){
    if(class(csceo$usage$descriptionText)=="character"){
      csceo$text <- csceo$usage$descriptionText
    } else {
      warning("'descriptionText' is not a character vector")
    }
  }
  
  
  # *Attaching CELLS2KEEP to output ----
  # If the SCE object contains over 50,000 cells, a random sample of 50,000 cells will be chosen for visualization in the application. 
  # Please note that all calculations are already completed, and this step is solely for promoting smooth and efficient visualization.
  
  if(ncol(csceo$SCE)>csceo$usage$nSubCells){
    if(verbose) message('Creating cell subset for plots because it exceeds ', csceo$usage$nSubCells, ' cells')
    csceo$CELLS2KEEP <- subsampling_func(csceo$SCE, cellsToKeep = csceo$usage$cells2keep, nmaxcell = csceo$usage$nSubCells)
    if(verbose) message('Finished')
  } else {
    csceo$CELLS2KEEP <- "all"
  }
  
  message('You can now launch the shiny server:\n e.g. launch_scX(cseo, point.size = 50, port = 9099, host = "0.0.0.0", launch.browser = F)')
  ##Free memory
  gc() 
  return(csceo)
}
