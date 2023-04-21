#### General ----

#### Upload a list of genes ----
#' @keywords internal
#' @noRd
genesList <- function(dataPath){
  if(file_ext(dataPath$name) == "txt"){
    #lista de genes txt separado por comas o newline
    genes <- readr::read_csv(file=dataPath$datapath, col_names = FALSE, col_types = readr::cols())
  }
  else{
    #lista de genes xlsx
    genes <- readxl::read_excel(path = dataPath$datapath, col_names = FALSE)
  }
  genes <- as.vector(as.matrix(genes))
  
  #saco los espacios, separo en "\\(" y saco "\\)"
  genes <- sub("\\)", "", unlist(strsplit(sub(" ", "", genes), split = "\\(")))
  
  return(unique(genes))
}

##### Correlation Boxes ----
# Cajitas de Luz ----
#' @keywords internal
#' @noRd
make_box <- function(ssce, selected.cells){
  # Create a vector that have a selected if the cell were selected
  ncells <- ncol(ssce)
  cells.v <- rep("unselected", ncells);names(cells.v) <- colnames(ssce)
  cells.v[selected.cells] <- "selected"
  names(cells.v) <- NULL
  # add the vector to SCE as a factor
  ssce$cajitasdeluz <- cells.v
  ssce$cajitasdeluz <- as.factor(ssce$cajitasdeluz)
  clusterss <- "cajitasdeluz"
  lengths <- table(colData(ssce)[,clusterss])
  ies <- rep(0, sum(lengths))
  contador <- 1
  for (cluster in names(lengths)){
    idx <- which(colData(ssce)[,clusterss] == cluster)
    ies[contador:(contador+length(idx)-1)] <- idx
    contador <- contador+length(idx)
  }
  jotas <- unlist(mapply(function(x,y){rep(x, y)},
                         seq(length(lengths)), lengths)
           )
  paraCorr <- Matrix::sparseMatrix(i = ies,
                           j = jotas,
                           x = 1,
                           dims = c(ncol(ssce), length(levels(colData(ssce)[,clusterss])))
              )
  rownames(paraCorr) <- colnames(ssce)
  
  return(list(paraCorr,
              ssce))
}

#' @keywords internal
#' @noRd
generar_correlacion <- function(mtx,mtx.cor,ssce, cluster){
  correlacion <- qlcMatrix::corSparse(t(mtx), mtx.cor[colnames(mtx),])
  colnames(correlacion) <- names(table(colData(ssce)[,cluster]))
  rownames(correlacion) <- rownames(mtx)
  correlacion[is.na(correlacion)] <- 0
  return(correlacion)
}

#' @keywords internal
#' @noRd
cajitasdeluz <- function(ssce, selected.cells, corr = 0.7){
  if(!("logcounts" %in% names(assays(ssce)))){
    stop("No 'logcounts' in the sce object")
  }
  lista <- make_box(ssce, selected.cells)
  correlacionMarker <- generar_correlacion(mtx = logcounts(lista[[2]]),
                                           mtx.cor = lista[[1]],
                                           ssce = lista[[2]],
                                           cluster = "cajitasdeluz")
  correlacionMarker <- data.frame(correlacionMarker)
  df <- correlacionMarker[correlacionMarker$selected > corr,1,drop=F]
  df <- df[order(df$selected, decreasing = TRUE),,drop=F]
  names(df)[1] <- "box.cor"
  return(df)
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