# Aesthetic, plot and visual functions

# Cell order from a partition to make the SpikePlot ----
#' @import RColorBrewer
#' @noRd
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

# CoExpression Tab ----
# CoExpression Palette ----
#' @noRd
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

# GridPlot ----
#' @noRd
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

# Heatmaps Functions ----
# vertical and horizontal Line in heatmaps
#' @noRd
vline <- function(x = 0){
  list(type = "line",
       y0 = -0.02, y1 = 1.02, yref = "paper",
       x0 = x, x1 = x,
       line = list(color = 'gray',
                   width = 0.8,
                   opacity = 0.8,
                   dash = 'dot')
       )
}

#' @noRd
hline <- function(y = 0){
  list(type = "line",
       x0 = -0.02, x1 = 1.02, xref = "paper",
       y0 = y, y1 = y,
       line = list(color = 'gray',
                   width = 0.8,
                   opacity = 0.8,
                   dash = 'dot')
       )
}

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


# Partition Tab ----
#' @noRd
set_val <- function(tab){
  val_fun <- function(j, i, x, y, w, h, col) { # add text to each grid
    grid::grid.text(round(tab,digits = 2)[i, j], x, y,gp= gpar(col="white"))
  }
  return(val_fun)
}


# MultiPlots Tab ----
# By clusters, plot the points with true values on the top
#' @noRd
reducedDimPlot_cluster <- function(sce,reducedDim,partition,cluster,alpha=0.5,palette="red"){
  df <- data.frame(x=reducedDim(sce,reducedDim)[,1],
                   y=reducedDim(sce,reducedDim)[,2],
                   value=factor((colData(sce)[,partition] %in% cluster),levels =c(TRUE,FALSE))
  )
  levels(df$value) <- c(cluster,"Other")
  g <- ggplot(df %>% arrange(desc(value)),aes(x=x,y=y,col=value)) + geom_point(alpha=alpha)+ 
    theme_classic() +  
    #labs(color=partition) + xlab(reducedDim) + ylab(reducedDim) +
	xlab(reducedDim) + ylab(reducedDim) + ggtitle(paste("Colored:", cluster)) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=c(palette,"lightgrey"))
  return(g)
}