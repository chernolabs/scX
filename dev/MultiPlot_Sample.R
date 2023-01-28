
sce <-list.sce$d1

lista <- list()
for(i in c("Week1","Week2","Week4","Week8")){
sce$value <- factor(colData(sce)[,"SampleID"] %in% c(i),levels = c(TRUE,FALSE))
lista[[i]]<- plotReducedDim(object = sce, dimred = "TSNE2D",
               ncomponents = 2,
               colour_by= "value") + 
  scale_color_manual(values=c("red","lightgrey")) + 
  labs(color=i)
  }

i <- "Week1"

g <- plotReducedDim(object = sce, dimred = "TSNE2D",
               ncomponents = 2,
               colour_by= "value") + 
  scale_color_manual(values=c("red","lightgrey")) + 
  labs(color=i)

g$mapping <- 0.0001


sem_Vtor <- levels(sce$SampleID)
lista <- list()
for(i in sem_Vtor){
  lista[[i]] <- reducedDimPlot_cluster(sce = sce,reducedDim = "TSNE2D",partition = "SampleID",cluster = i,alpha = 0.5,palette = )
}
plot_grid(plotlist =lista)
rank(factor(colData(sce)[,"value"]),ties.method="first")
ggplot() + geom_point(aes(x=reducedDim(sce,"TSNE2D")[,1],y=reducedDim(sce,"TSNE2D")[,2],
                          col=factor(colData(sce)[,"value"])),alpha=0.1) + 
  theme_classic() +  
  scale_color_manual(values=c("red","lightgrey")) + 
  labs(color=i)
               