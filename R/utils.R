#' @import shiny
#' @import shinydashboard
#' @import shinyWidgets
#' @import shinycssloaders
#' @import shinycustomloader
#' @import shinydisconnect
#' @import shinyFeedback
#' @import shinyjs
#' @import shinyBS
#' @import DT
#' @import Matrix
#' @import reshape2
#' @import plotly
#' @import ggplot2
#' @import ComplexHeatmap
#' @import magick
#' @importFrom magrittr %>%
#' @importFrom htmltools HTML
#' @importFrom htmltools tags
#' @importFrom scales rescale
#' @import dplyr
#' @import tidyr

# General functions

# Upload a list of genes ----
#' @keywords internal
#' @noRd
genesList <- function(dataPath){
  filext <- tools::file_ext(dataPath$name)
  if(filext %in% c("txt","text","csv")){
    # txt list of genes, separated by commas or newline
    genes <- readr::read_csv(file=dataPath$datapath, col_names = FALSE, col_types = readr::cols())
  } else if(filext %in% c("xls","xlsx")){
    # xlsx list of genes
    genes <- readxl::read_excel(path = dataPath$datapath, col_names = FALSE)
  } else {
    genes <- NULL
	warning('Supported file formats for gene lists are .txt, .csv, .xls, .xlsx')
  }
  genes <- as.vector(as.matrix(genes))
  
  # remove spaces, separate at "\\(" and remove "\\)"
  genes <- sub("\\)", "", unlist(strsplit(sub(" ", "", genes), split = "\\(")))
  
  return(unique(genes))
}

# Correlation Boxes ----
# Cajitas de Luz ----
#' @keywords internal
#' @noRd
make_box <- function(ssce, selected.cells){
  # Create a vector that has 'selected' for selected cells
  ncells <- ncol(ssce)
  cells.v <- rep("unselected", ncells);names(cells.v) <- colnames(ssce)
  cells.v[selected.cells] <- "selected"
  names(cells.v) <- NULL
  # add the vector to SCE as a factor
  ssce$selectionbox <- cells.v
  ssce$selectionbox <- as.factor(ssce$selectionbox)
  
  # Note: this code was made for partitions with more than 2 levels
  # For just 2 levels (selected/unselected), it's redundant but still works
  partition <- "selectionbox"
  lengths <- table(colData(ssce)[,partition])
  ies <- rep(0, sum(lengths))
  contador <- 1
  for (cluster in names(lengths)){
    idx <- which(colData(ssce)[,partition] == cluster)
    ies[contador:(contador+length(idx)-1)] <- idx
    contador <- contador+length(idx)
  }
  jotas <- unlist(mapply(function(x,y){rep(x, y)},
                         seq(length(lengths)), lengths)
           )
  paraCorr <- Matrix::sparseMatrix(i = ies,
                           j = jotas,
                           x = 1,
                           dims = c(ncells, length(lengths)),
						   dimnames = list(colnames(ssce), names(lengths))
              )
			  
  return(list(boxes=paraCorr,
              sce=ssce))
}

#' @keywords internal
#' @noRd
box_correlation <- function(ssce, selected.cells, corr = 0.7){
  if(!("logcounts" %in% assayNames(ssce))){
    stop("No 'logcounts' in the sce object")
  }
  lista <- make_box(ssce, selected.cells)
  
  correlationMarker <- sparseCor(t(logcounts(lista[["sce"]])),
								lista[["boxes"]])
  correlationMarker[is.na(correlationMarker)] <- 0
  
  correlationMarker <- data.frame(correlationMarker)
  df <- correlationMarker[correlationMarker$selected > corr,1,drop=F]
  df <- df[order(df$selected, decreasing = TRUE),,drop=F]
  names(df)[1] <- "box.cor"
  return(df)
}

#' @keywords internal
#' @noRd
sparseColMax <- function(X){
  X <- as(X,"dgCMatrix")
  sapply(seq(ncol(X)), function(k){
    p0 <- X@p[k]
    pf <- X@p[k+1]
    if(p0==pf){ # column k has no elements
      0
    }else{
      max(X@x[(p0+1):pf])
    }
  })
  # as(*, "sparseVector") unnecessary?
}

#' @keywords internal
#' @noRd
sparseRowMax <- function(X){
  X <- t(X)
  sparseColMax(X)
}

#' @keywords internal
#' @noRd
sparseCor <- function(X,Y=NULL){
  X <- as(X,"dgCMatrix")
  N <- nrow(X)
  
  Xmeans <- colMeans(X)
  
  if(!is.null(Y)){
    if(nrow(Y)!=N) stop("X and Y have different number of rows")
    Y <- as(Y,"dgCMatrix")
    Ymeans <- colMeans(Y)
    
    # almost standard deviation, (n-1) missing
    Xsd <- sqrt(colSums(X^2) - N*Xmeans^2)
    Ysd <- sqrt(colSums(Y^2) - N*Ymeans^2)
    
    # almost cov, (n-1) missing, would cancel out with Xsd*Ysd:
    covn1 <- as.matrix(crossprod(X, Y)) - N*tcrossprod(Xmeans, Ymeans) # # t(X) %*% Y - N * (Xmeans %o% Ymeans)
    return(covn1/tcrossprod(Xsd, Ysd)) # (Xsd %o% Ysd)
  } else {
    # almost cov, (n-1) missing 
    covn1 <- as.matrix(crossprod(X)) - N*tcrossprod(Xmeans) # t(X) %*% X - N * (Xmeans %o% Xmeans) 
    Xsd <- sqrt(diag(covn1)) # almost standard dev, (n-1) missing, would cancel out with cov
    return(covn1/tcrossprod(Xsd)) # (Xsd %o% Xsd)
  }
}

#CoExpression DF ----
#' @keywords internal
#' @noRd
COexp_Vtor <- function(sce,genes){
  coexp <- assay(sce,"logcounts")[genes,] %>% apply(1,function(x){ifelse(x>0,"1","0")})
  vtor  <- setNames(c("None",genes[2],genes[1],"Both"),
                    nm = c("00","01","10","11"))
  vt    <- factor(vtor[paste0(coexp[,1],coexp[,2])],
                  levels = c("None",genes[1],genes[2],"Both"))
  return(vt)
}

### Add imagenes to shiny
#' @keywords internal
#' @noRd
.onAttach <- function(libname, pkgname) {
  shiny::addResourcePath('www',
                         system.file('www',
                                     package = 'scX'))
}

#' @keywords internal
#' @noRd
toWebGL <- function(p) {
	# 06/10/2023 modified from:
	# https://github.com/zeehio/plotly.R/commit/d9652a6f55d6969037a18983dd10d59a8b225886#diff-c0e4d1c45c2e1d458c837daac8876c3071b5c4794c29fcba01e3b250a8932df8
	# check for plotly updates / merges
    if (ggplot2::is.ggplot(p)) {
        p <- plotly::plotly_build(p)
	}
	traces_without_hoveron <- plotly:::glTypes()
	trace_idx <- vapply(
		p$x$data,
		function(trace) trace$type %in% traces_without_hoveron,
		logical(1)
	)
	p <- style(p, hoveron = NULL, traces = which(trace_idx))
    p$x$.plotlyWebGl <- TRUE
    p
}