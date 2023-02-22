#### General ----
###Cell order from a partition to make the SpikePlot ----
Col.and.Order <- function(partition,sce){
  #Color
  lvs <- levels(colData(sce)[,partition])
  
  cellTypeCol <- c(brewer.pal(12, "Set3"),
                   brewer.pal(8, "Dark2"),
                   brewer.pal(8, "Accent"),
                   brewer.pal(9, "Pastel1"))[3:(2+length(lvs))]
  names(cellTypeCol) <- lvs
  
  #Order
  ord <-order(ordered(colData(sce)[,partition]))
  
  return(list(colPart = cellTypeCol,
              ordPart =ord))
}

#### Upload a list of genes ----

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

##### Correlation Boxs ----

require(DropletUtils)
require(SingleCellExperiment)
require(scater)
require(scran)
require(Matrix)

# Cajitas de Luz ----

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
  paraCorr <- sparseMatrix(i = ies,
                           j = jotas,
                           x = 1,
                           dims = c(ncol(ssce), length(levels(colData(ssce)[,clusterss])))
              )
  rownames(paraCorr) <- colnames(ssce)
  
  return(list(paraCorr,
              ssce))
}

generar_correlacion <- function(mtx,mtx.cor,ssce, cluster){
  require(qlcMatrix)
  correlacion <- corSparse(t(mtx), mtx.cor[colnames(mtx),])
  colnames(correlacion) <- names(table(colData(ssce)[,cluster]))
  rownames(correlacion) <- rownames(mtx)
  correlacion[is.na(correlacion)] <- 0
  return(correlacion)
}

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

#### Markers Tab ----
####  Scientifict notation DT ----
# js <- c(
#   "function(row, data, displayNum, index){",
#   "  var x = data[3];",
#   "  $('td:eq(3)', row).html(x.toExponential(2));",
#   "  var y = data[4];",
#   "  $('td:eq(4)', row).html(y.toExponential(2));",
#   "}"
# )

#### Diff Expression Tab ----
####  Scientifict notation DT ----
js_volcano <- c(
  "function(row, data, displayNum, index){",
  "  var x = data[2];",
  "  $('td:eq(2)', row).html(x.toExponential(2));",
  "  var y = data[3];",
  "  $('td:eq(3)', row).html(y.toExponential(2));",
  "}"
)


### CoExpression Tab----
##CoExpression Palette  ----
CoExp_Col  <- function(genes,sce){
  foo <- assay(sce,"logcounts")[genes,] %>%  
    apply(1,function(x){ #Norm by gene
      (x- mean(x)) / sd(x)
    }) %>% apply(2,function(x){ #Scale
      rescale(x)
    })
  
  foo1 <- as.numeric(foo[,1])
  foo2 <- as.numeric(foo[,2])
  radio <- sqrt(foo1^2 + foo2^2)
  angle <- atan(foo2/foo1)
  angle[which(is.nan(angle))] <- 0
  col <- hsv(h = 1-2*angle/(3*pi),
             s = (mapply(min, 1, radio)+0.1)/1.1,
             v = (mapply(min, 1, radio)+0.5)/1.5,
             alpha = (mapply(max, foo1, foo2)+1)/2
         )
  return(col)
}

#GridPlot ----

plot2dgradient <- function(res = 50, gen1 = "Gen1", gen2 = "Gen2"){
  par(mfrow=c(1,1),
      mar=c(2,2,1,1),
      mgp=c(1,0,0))
  
  grid <- c(0, seq(res)/res) # Square limits
  # left/bottom is grid except last element;
  #right/top is grid except first element;
  
  colseq <- c(0, seq(res-1)/(res-1)) 
  
  col1 <- rep(colseq, res)
  col2 <- sort(rep(colseq, res))
  alpha <- (mapply(max, col1, col2)+1)/2 # Scale alpha 0.5 to 1
  radio <- sqrt(col1^2 + col2^2)
  angulo <- atan(col2/col1)
  angulo[which(is.nan(angulo))] <- 0
  hue <- 1-2*angulo/(3*pi)
  
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = F, xlab = gen1, ylab = gen2)
  rect(xleft = rep(grid[1:length(grid)-1], res),
       ybottom = sort(rep(grid[1:length(grid)-1], res)),
       xright = rep(grid[2:length(grid)], res),
       ytop = sort(rep(grid[2:length(grid)], res)),
       border = NA,
       col = hsv(h = hue,
                 s = (mapply(min, 1, radio)+0.1)/1.1,
                 v = (mapply(min, 1, radio)+0.5)/1.5,
                 alpha = alpha
             )
  )
  text(0, 1, '1', pos = 2)
  text(0, 0, '0', pos = 2)
  text(1, 0, '1', pos = 4)
}

#CoExpression DF ----
COexp_Vtor <- function(sce,genes){
  coexp <- assay(sce,"logcounts")[genes,] %>% apply(1,function(x){ifelse(x>0,"1","0")})
  vtor  <- setNames(c("None",genes[2],genes[1],"Both"),
                    nm = c("00","01","10","11"))
  vt    <- factor(vtor[paste0(coexp[,1],coexp[,2])],
                  levels = c("None",genes[1],genes[2],"Both"))
  return(vt)
}


#### Heatmaps Functions -----
# 
# #Line in heatmaps
# vline <- function(x = 0){
#   list(type = "line", 
#        y0 = -0.02, y1 = 1.02, yref = "paper",
#        x0 = x, x1 = x, 
#        line = list(color = 'gray',
#                    width = 0.8,
#                    opacity = 0.8,
#                    dash = 'dot')
#        )
# }
# 
# hline <- function(y = 0){
#   list(type = "line", 
#        x0 = -0.02, x1 = 1.02, xref = "paper",
#        y0 = y, y1 = y, 
#        line = list(color = 'gray',
#                    width = 0.8,
#                    opacity = 0.8,
#                    dash = 'dot')
#        )
# }
# 
# 
# ParamHeatm <- function(partition,sce){
#   Ncells   <- table(colData(sce)[,partition])
#   lineVals <- cumsum(Ncells)[1:(length(Ncells)-1)]
#   tickVals <- Ncells/2
#   for (i in 2:length(tickVals)){
#     tickVals[i] <- tickVals[i] + lineVals[i-1]
#   }
#   
#   foo <- list(tickVals = tickVals,
#               lineVals=lineVals,
#               lev = names(Ncells))
#   return(foo)
# }


##### Partition Tab ----

set_val <- function(tab){
  val_fun <- function(j, i, x, y, w, h, col) { # add text to each grid
    grid.text(round(tab,digits = 2)[i, j], x, y,gp= gpar(col="white"))
  }
  return(val_fun)
}




##### MultiPlots ----
#By clusters, plot the points with true values on the top
reducedDimPlot_cluster <- function(sce,reducedDim,partition,cluster,alpha=0.5,palette="red"){
  df <- data.frame(x=reducedDim(sce,reducedDim)[,1],
                   y=reducedDim(sce,reducedDim)[,2],
                   value=factor((colData(sce)[,partition] %in% cluster),levels =c(TRUE,FALSE))
  )
  levels(df$value) <- c(cluster,"Other")
  g <- ggplot(df %>% arrange(desc(value)),aes(x=x,y=y,col=value)) + geom_point(alpha=alpha)+ 
    theme_classic() +  
    labs(color=partition) + xlab(reducedDim) + ylab(reducedDim) +
    scale_color_manual(values=c(palette,"lightgrey"))
  return(g)
}

### SCE function ----
ldf_func <- function(sce, partition){
  
  cat(partition, ":\n")
  # le saco los factores vacios
  colData(sce)[,partition] <- droplevels(colData(sce)[,partition])
  
  # calculo lfmrk ----
  type           <- c('all','any','some')
  lfmrk <- list()
  for(itype in seq_along(type)){
    lfmrk[[type[itype]]]    <- findMarkers(sce,groups=colData(sce)[,partition],full.stats=TRUE,
                                           row.data=rowData(sce)[,c(2,9)],
                                           direction='up',pval.type=type[itype])
  }
  
  # calculo ldf ----
  ldf_t<-list()
  numCores <- 30
  minSize  <- 20
  registerDoParallel(numCores)
  
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
    
    # (6.1.1) Calculo de entropia ----
    # . comienza el calculo de entropia 
    # para la cuenta de entropia solo analizo clusters de mas de minSize celulas
    lab   <- colData(sce)[,partition]
    tt    <- table(lab)
    labOK <- tt > minSize
    labOK <- names(labOK)[labOK]
    
    Z   <- logcounts(sce)[a4,]  
    
    pattern <- rep(0,ncol(Z))
    pattern[colData(sce)[,partition]%in%coi] <- 1
    boxcor<-apply(Z,1,function(x){
      return(cor(x,pattern))
    })
    
    #calcula entropio/informatividad de cada gen
    resS <- c()
    #for(igen in 1:nrow(Z)){
    aa<-foreach(igen=1:nrow(Z),.combine='c') %dopar% {
      a     <- aggregate(Z[igen,lab%in%labOK],by=list(lab[lab%in%labOK]),mean)
      #x <- (a[,2] - min(a[,2]))/(sum(a[,2])-min(a[,2])*nrow(a))
      x <- a[,2]/sum(a[,2])
      s <- -sum(x[x>0]*log2(x[x>0]))
      s <- s/log2(length(x))
      s.max <- a[which.max(x),1]
      
      muCI  <- a[a[,1]==coi,2]
      muCI2 <- sort(a[,2],decreasing=TRUE)[2]
      resS <- c(resS,c(1-s,s.max,muCI,muCI2))
    }
    
    mInfo <- data.frame(matrix(aa,byrow=TRUE,ncol=4))
    colnames(mInfo) <- c("selectivity","clusterMax",'meanX','meanX2nd')
    rownames(mInfo) <- rownames(Z)
    mInfo[,2] <- names(lfmrk[[1]])[mInfo[,2]]
    
    
    
    df <- data.frame(boxcor=boxcor,mInfo[,c("selectivity", "meanX", "clusterMax", "meanX2nd")],robustness=a[names(mrkrs)])
    
    #me quedo con los marcadores que tiene pico de expresion en el lcuster de interes
    df <- df[df$clusterMax==names(lfmrk[[1]])[ic],]
    df <- df[order(df$boxcor,decreasing=TRUE),]
    
    ldf_t[[coi]]<-df
  }
  return(ldf_t)
}



