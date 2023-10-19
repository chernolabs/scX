# Aesthetic, plot and visual functions

# Cell order from a partition to make the SpikePlot ----
#' @import RColorBrewer
#' @keywords internal
#' @noRd
Col.and.Order <- function(partition,sce){
  #Color
  lvs <- levels(colData(sce)[,partition])
  L <- length(lvs)
  
  cellTypeCol <- c(brewer.pal(12, "Set3"),
                   brewer.pal(8, "Dark2"),
                   brewer.pal(8, "Accent"),
                   brewer.pal(9, "Pastel1"))
	if (L<length(cellTypeCol)-2){
		cellTypeCol <- cellTypeCol[3:(2+L)]
	}else{
		cellTypeCol <- rainbow(L)
	}
  names(cellTypeCol) <- lvs
  
  #Order
  ord <-order(ordered(colData(sce)[,partition]))
  
  return(list(colPart = cellTypeCol,
              ordPart =ord))
}

# Legend columns number and width for SpikePlot horizontal legend ----
#' @param names Vector of legend names
#' @param L Numeric, width of plot (i.e. xlim or maximum x value)
#' @keywords internal
#' @noRd
legend_col <- function(names, L){
  # `legend` allows ncol parameter (not nrow) but fills spaces by rows (each row in col1, then col2, etc)
  # hence, fixing ncol might result in empty columns. this function optimizes ncol to avoid this.
  N <- length(names)
  nameL <- strwidth(names)
  
  nrows <- 1
  ncols <- N
  rem <- 0
  textwidth <- nameL
  totalwidth <- sum(textwidth)
  while(totalwidth>L && nrows<N){ #if the sum of column widths is larger than plot width
    nrows <- nrows+1 # a rows is added
    rem <- N%%nrows
    ncols <- ifelse(rem==0, N%/%nrows, N%/%nrows+1) # no. of columns depends on whether there's a remainder
    width_mat <- matrix(c(nameL, rep(0, ncols*nrows-N)), ncol=ncols, nrow=nrows) # complete matrix with 0s if necessary
    textwidth <- apply(width_mat, 2, max) #width for each column
    totalwidth <- sum(textwidth)
  }
  return(list(ncol=ncols, colwidth=textwidth))
}


# CoExpression Tab ----
# CoExpression Palette
#' @keywords internal
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

# GridPlot
#' @keywords internal
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
#' @keywords internal
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

#' @keywords internal
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
#' @keywords internal
#' @noRd
set_val <- function(tab){
  val_fun <- function(j, i, x, y, w, h, col) { # add text to each grid
    grid::grid.text(round(tab,digits = 2)[i, j], x, y,gp= grid::gpar(col="white"))
  }
  return(val_fun)
}


# MultiPlots Tab ----
# By clusters, plot the points with true values on the top
#' @keywords internal
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

#### Fields Tab ----

#Modified from https://rdrr.io/github/davismcc/scater/src/R/plotDots.R to allow column vectors
plotDots_fields <- function(df, features, group = NULL, block=NULL,
                            exprs_values = "logcounts", detection_limit = 0, zlim = NULL, 
                            colour = color, color = NULL,
                            max_detected = NULL, other_fields = list(),
                            by_exprs_values = exprs_values,
                            swap_rownames = NULL,
                            center = FALSE,
                            scale = FALSE,
                            assay_name=exprs_values,
                            by_assay_name=by_exprs_values)
{
  
  # if (is.null(group)) {
  #   group <- rep("all", ncol(object))
  # } else {
  #   group <- scater::retrieveCellInfo(object, group, search="colData")$value
  #   noAll <- T
  # }
  # object <- .swap_rownames(object, swap_rownames)
  # features <- .handle_features(features, object)
  summarized1 <- df %>% group_by(across(all_of(group))) %>% summarise(ave = across(all_of(features),mean,na.rm = TRUE),
                                                                      prop = across(all_of(features),function(x){sum(x>0,na.rm = T)/length(x)})
  )
  
  ave <- summarized1$ave %>% t
  ave[is.na(ave)] <- 0
  num <- summarized1$prop %>% t
  ave[is.na(num)] <- 0
  group.names <- if(!is.null(group)){summarized1[,group,drop=T]} else{"All"}
  
  #Function from https://rdrr.io/github/davismcc/scater/src/R/plotHeatmap.R#sym-.heatmap_scale
  .heatmap_scale <- function(x, center, scale, colour=NULL, zlim=NULL, symmetric=NULL) {
    
    if (center & !is.null(group)) {
      x <- x - rowMeans(x)
    }
    if (scale & !is.null(group)) {
      if (!center & any(rowSums(x) == 0)) {
        stop("Cannot include non-expressed genes when scale=TRUE.")
      }
      x <- x / sqrt(rowSums(x^2) / (ncol(x) - 1))
    }
    if (is.null(zlim)) {
      if (center & !is.null(group)) {
        extreme <- max(abs(x))
        zlim <- c(-extreme, extreme)
      } else {
        zlim <- range(x)
      }
    }
    if (is.null(colour)) {
      if (center) {
        colour <- rev(RColorBrewer::brewer.pal(9, "RdYlBu"))
      } else {
        colour <- viridis::viridis(9)
      }
    }
    x[x < zlim[1]] <- zlim[1]
    x[x > zlim[2]] <- zlim[2]
    list(
      x = x,
      colour = colour,
      colour_breaks = seq(zlim[1], zlim[2], length.out=length(colour) + 1L),
      colour_scale = scale_colour_gradientn(colours = colour, limits = zlim),
      zlim = zlim
    )
  }
  heatmap_scale <- .heatmap_scale(ave, center=center, scale=scale, colour=colour, zlim=zlim)
  
  # Creating a long-form table.
  evals_long <- data.frame(
    Feature=rep(features, ncol(num)),
    Group=rep(group.names, each=nrow(num)),
    NumDetected=as.numeric(num),
    Average=as.numeric(heatmap_scale$x)
  )
  # if (!is.null(max_detected)) {
  #   evals_long$NumDetected <- pmin(max_detected, evals_long$NumDetected)
  # }
  
  # Adding other fields, if requested.
  # vis_out <- .incorporate_common_vis_row(evals_long, se = object,
  #                                        colour_by = NULL, shape_by = NULL, size_by = NULL, 
  #                                        by_assay_name = by_assay_name, other_fields = other_fields,
  #                                        multiplier = rep(.subset2index(features, object), ncol(num)))
  # evals_long <- vis_out$df
  ggplot(evals_long) + 
    geom_point(aes(x=.data$Group, y=.data$Feature, size=.data$NumDetected, col=.data$Average)) +
    scale_size(limits=c(0, max(evals_long$NumDetected))) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(size=0.5,colour = "grey80"),
      panel.grid.minor = element_line(size=0.25,colour = "grey80")) +
    heatmap_scale$colour_scale 
}