server <- function(input, output,session) {
  observeEvent(input$button.config,{
    QC_Server("qc",sce = cseo$SCE)  
    markersServer(id="markers",sce=cseo$SCE,ldf = cseo$ldf,point.size = input$point.size)
    N_markersServer(id="n_markers",sce=cseo$SCE,point.size = input$point.size)
    ExpressionServer(id="Exp",sce=cseo$SCE,point.size = input$point.size)
    COExpServer(id="Co-exp",sce=cseo$SCE,point.size = input$point.size)
    VolcanoServer(id="volcano",sce=cseo$SCE,sce.markers = cseo$sce.markers)
    VT_Server(id = "tools",sce =cseo$SCE)
    MultiPlotsServer(id = "MP",sce =cseo$SCE)
    Clusters_Server("cluster",sce = cseo$SCE)                                                                                                                                                                                                     
  })
}
