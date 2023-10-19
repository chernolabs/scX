#Maeke correlation betw
col1 <- "SampleID"
val<-c("nCounts","nFeatures")
df <- colData(sce) %>% as.data.frame() %>% group_by(SampleID) %>% summarise(Correlation = cor(cbind(across(all_of(val)))))
rownames(df$Correlation)  <- paste0(df[,col1,drop=T],"_",rownames(df$Correlation))
colnames(df$Correlation) <- val

Heatmap(t(df$Correlation),
        column_split = ,
        show_row_names = F)

dta<-t(df$Correlation)
Heatmap(dta,
        col = if(max(dta)==min(dta)) {viridis(1)} else {viridis(100)},
        border =F,
        name = "Gene expression",
        cluster_rows =T,
        cluster_columns = T,
        row_names_side = "left",
        column_names_rot = 45,
        # column_title_rot = 45,
        # column_title_gp = gpar(fontsize = 10),
        column_title_side = "bottom",
        # column_title = "Partition",
        row_title = "Genes",
        row_gap = unit(1, "mm"),
        column_gap = unit(1, "mm"),
        show_row_names = T,
        # show_column_names = T,
        # top_annotation = ht,
        column_split = df[,col1,drop=T],
        column_labels = sapply(str_split(colnames(dta),pattern = "_"),function(x){tail(x,n=1)}),
        cluster_column_slices = F
        # use_raster = TRUE,
        # raster_by_magick = TRUE
)
if(input$partitionColor != 'None'){
        ht <-   HeatmapAnnotation(Type = df()[,input$partitionColor],
                                  col=list(Type=OrderPartReact()$colPart),
                                  show_legend = c(Type =FALSE),
                                  annotation_label = c(input$partitionColor),
                                  show_annotation_name = T)
        col_split <- df()[,input$partitionColor]
        } else{
          ht <- NULL
          col_split <- NULL
        }
Heatmap(t(df$Correlation),
        col = if(max(t(df$Correlation))==min(t(df$Correlation))) {viridis(1)} else {viridis(100)},
        border =F,
        name = "Gene expression",
        cluster_rows = T,
        cluster_columns = T,
        row_names_side = "left",
        # column_names_rot = 45,
        column_title_rot = 45,
        # column_title_gp = gpar(fontsize = 10),
        column_title_side = "bottom",
        # column_title = "Partition",
        row_title = "Genes",
        row_gap = unit(1, "mm"),
        column_gap = unit(1, "mm"),
        show_row_names = T,
        show_column_names = T,
        # top_annotation = ht,
        column_split = df[,col1,drop=T],
        cluster_column_slices = F,
        use_raster = F,
        raster_by_magick = F
)

col1 <- "None"
df <- colData(sce) %>% as.data.frame() %>% summarise(Correlation = cor(cbind(across(all_of(val))))) %>% .[["Correlation"]]

###Dotplot to columns ----

#Modified from https://rdrr.io/github/davismcc/scater/src/R/plotDots.R to allow column vectors
plotDots_fields <- function(object, features, group = NULL, block=NULL,
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
  
  if (is.null(group)) {
    group <- rep("all", ncol(object))
  } else {
    group <- retrieveCellInfo(object, group, search="colData")$value
  }
  
  # object <- .swap_rownames(object, swap_rownames)
  # features <- .handle_features(features, object)
  group <- factor(group)
  
  # Computing, possibly also batch correcting.
  ids <- DataFrame(group=group)
  if (!is.null(block)) {
    ids$block <- retrieveCellInfo(object, block, search="colData")$value
  }
  
  summarized <- summarizeAssayByGroup(
    as.matrix(t(colData(object)[,as.character(features), drop = FALSE])),
    ids=ids, statistics=c("mean", "prop.detected"),
    threshold=detection_limit)
  
  ave <- assay(summarized, "mean")
  num <- assay(summarized, "prop.detected")
  group.names <- summarized$group
  
  if (!is.null(block)) {
    ave <- correctGroupSummary(ave, group=summarized$group, block=summarized$block)
    num <- correctGroupSummary(num, group=summarized$group, block=summarized$block, transform="logit")
    group.names <- factor(colnames(ave), levels = levels(summarized$group))
  }
  #Function from https://rdrr.io/github/davismcc/scater/src/R/plotHeatmap.R#sym-.heatmap_scale
  .heatmap_scale <- function(x, center, scale, colour=NULL, zlim=NULL, symmetric=NULL) {
    
    if (center) {
      x <- x - rowMeans(x)
    }
    if (scale) {
      if (!center & any(rowSums(x) == 0)) {
        stop("Cannot include non-expressed genes when scale=TRUE.")
      }
      x <- x / sqrt(rowSums(x^2) / (ncol(x) - 1))
    }
    if (is.null(zlim)) {
      if (center) {
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
  if (!is.null(max_detected)) {
    evals_long$NumDetected <- pmin(max_detected, evals_long$NumDetected)
  }
  
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

plotDots(object = sce,
         features = c("nCounts","nFeatures"),
         group ="ds1clusterP",scale = T,center = T)



### COexp-----
genes <- c(sample(rownames(sce),10),"Sema3a","Nell1")
test <- logcounts(sce)[genes,]
test <- test %>% apply(1,function(x){x>0})
# df_t <- combn(x = test[2,],m = 2,FUN = function(x){sort(names(x))},simplify = F)
mm <- apply(test,1,function(w){
  df_t1 <- combn(x = w,m = 2,simplify = F)
  df <-lapply(df_t1,function(x){
    x <- x[sort(names(x))]
    c(names(x),paste0(as.numeric(x),collapse = ""))
  }) %>% do.call("rbind",.) %>% as.data.frame
})
fnl <- do.call("rbind",mm) %>% group_by(V1,V2,V3) %>% summarise(N=n())
fnl$V3 <- as.factor(fnl$V3)
levels(fnl$V3) <- c("None","01","10","Both")
fnl$V3 <- as.character(fnl$V3)
fnl[fnl$V3 == "01","V3"] <- fnl[fnl$V3 == "01","V2"]
fnl[fnl$V3 == "10","V3"] <- fnl[fnl$V3 == "10","V1"]

fnl %>% group_by(V1,V2)
