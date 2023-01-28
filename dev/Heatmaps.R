sce <-list.sce[[1]]

library(circlize)
if(T){
  
  col_pal <- c(brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))
  col_fun2 = colorRampPalette(colors = c("grey","blue"))
  df <- colData(sce) %>% as.data.frame
  
  Region_Cluster <- df[,"clustersOK"] %>% as.character
  Region_col      <- setNames(nm =unique(Region_Cluster),
                              object = col_pal[1:length(unique(Region_Cluster))])
  
  Sample_Cluster <- df[,"SampleID"]
  Sample_col      <- setNames(nm =unique(df[,"SampleID"]),object = col_pal[1:length(unique(df[,"SampleID"]))])
  
  LR_Cluster <- df[,"LR"]
  LR_col      <- setNames(nm =unique(LR_Cluster),object = col_pal[1:length(unique(LR_Cluster))])
  
  # Sex_Cluster <- metadata$Sex
  # Sex_Cluster[is.na(Sex_Cluster)]  <- 'NA'
  # Sex_col      <- setNames(nm =unique(Sex_Cluster),object = col_pal[1:length(unique(metadata$Sex))])
  # 
  # Age_Cluster <- metadata$Age
  # Age_col      <- setNames(nm = sort(unique(metadata$Age)),object = col_fun2(length(sort(unique(metadata$Age)))))
  # 
  # GF_Cluster <- metadata$GrowthFactorReceptor
  # GF_col      <- setNames(nm =unique(GF_Cluster),object = col_pal[1:length(unique(GF_Cluster))])
  # 
  
  
  
  hb <-  HeatmapAnnotation(df = df[,c("clustersOK","SampleID","LR")],
                       col=list(clustersOK=Region_col,
                                SampleID = Sample_col,
                                LR = LR_col))
mk <-list.ldf[[1]]$SampleID$new_8w %>% rownames %>% .[1:20]
} #HB (row) anottation
a <-Col.and.Order(partition = "SampleID",sce = sce)
col_fun = colorRamp2(c(0, 2), c("lightgrey","red"))
mat <- logcounts(sce)[mk,]  %>% as.matrix

orde <- T
Heatmap(mat, name = "mat", 
        col = col_fun,
        show_row_names = ifelse(nrow(mat)>20,F,T),
        show_column_names = F,
        top_annotation = hb,
        column_order = if(orde) {a$ordPart} else{NULL},
        cluster_columns = !orde)

