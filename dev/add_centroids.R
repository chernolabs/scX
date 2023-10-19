centroide <- F
size.l <- ncol(sce) < 200000
p <- plot_ly(type = "scatter3d", mode = "markers")
if(size.l){
  p <- p %>% add_trace(x = ~reducedDim(sce,"TSNE")[,1], y=~reducedDim(sce,"TSNE")[,2],
                       z=~reducedDim(sce,"TSNE")[,3],
                       text= ~colData(sce)[,"SampleID"],
                       hoverinfo = 'text',
                       color = ~logcounts(sce)["Nell1",],
                       name = ~colData(sce)[,"SampleID"],
                       size = I(20),span=I(0))
}
if(!size.l | centroide){
  p <- p %>% add_trace(x = ~tapply(reducedDim(sce,"TSNE")[,1],colData(sce)[,"SampleID"],mean), y=~tapply(reducedDim(sce,"TSNE")[,2],colData(sce)[,"SampleID"],mean),
                       z=~tapply(reducedDim(sce,"TSNE")[,3],colData(sce)[,"SampleID"],mean),
                       text= ~levels(colData(sce)[,"SampleID"]),
                       hoverinfo = 'text',
                       color = ~tapply(logcounts(sce)["Nell1",],colData(sce)[,"SampleID"],mean),
                       name = if(!size.l){~levels(colData(sce)[,"SampleID"])} else{"Centroid"},
                       size = I(20*5),span=I(0))
}
p %>% layout(dragmode = "select",
             scene = list(xaxis = list(title = 'Dim1',showgrid=F,visible=F),
                          yaxis = list(title = 'Dim2',showgrid=F,visible=F),
                          zaxis = list(title = 'Dim3',showgrid=F,visible=F)),
             legend= list(x=1,y=1),
             showlegend = TRUE,
             # title = ifelse(length(ExpressionF()$Genes)>1,'Mean Expression', 'Expression'),
             margin = list(l = 0,
                           r = 10,
                           b = 0,
                           t = 40,
                           pad = 0)) %>% 
  colorbar(title = "log(counts)",x=0,y=1) %>% 
  toWebGL()
