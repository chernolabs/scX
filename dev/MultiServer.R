server <- function(input, output,session) {
  
  #MultiDataset ----
  all <- reactive({
           list(sce=list.sce[[input$Dataset]],
                sce.markers = list.sce.markers[[input$Dataset]],
                ldf = list.ldf[[input$Dataset]]
            )
         })
  
  observeEvent(c(all()),ignoreInit = T,{
    sendSweetAlert(
      session = session,
      title = "Changing dataset",
      text = "This could take time",
      type = "info"
    )
    session$reload()
  })
  
  ### GlobalConfig ----
  observeEvent(c(all(),input$button.config),{
    #Summary Tab
    QC_Server("qc",sce = all()$sce)  
    
    #Markers Tab
    markersServer(id="markers",sce=all()$sce,ldf = all()$ldf,point.size = input$point.size)
    N_markersServer(id="n_markers",sce=all()$sce,point.size = input$point.size)
    
    #Gene Expression Tab
    ExpressionServer(id="Exp",sce=all()$sce,point.size = input$point.size)
    COExpServer(id="Co-exp",sce=all()$sce,point.size = input$point.size)
    
    #Diff Expression Tab
    VolcanoServer(id="volcano",sce=all()$sce,sce.markers = all()$sce.markers)
    
    #Partition Tab
    Clusters_Server("cluster",sce = all()$sce)
    
    #Tools Tab 
    VT_Server(id = "tools",sce =all()$sce)
    MultiPlotsServer(id = "MP",sce =all()$sce)
    
  })
}

