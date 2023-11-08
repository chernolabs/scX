#####                              ######## 
#####        Markers Tab           ######## 
#####                              ######## 

##### Marker UI Module ----
markersUI <- function(id = "markers") {
  tagList(
    useShinyjs(),
    fluidRow(
      column(4,
        box(title = htmltools::span(icon("gears"), " Settings"),
            width = NULL, status = "primary",solidHeader = T,collapsible = T,
          fluidRow(
            column(3,style='padding-left:12px; padding-right:3px;',align="center",
              pickerInput(NS(id,"DimType"), "  # dims", choices = NULL,width = NULL)
            ),
			column(4,style='padding-left:3px; padding-right:3px;',align="center",
              pickerInput(NS(id,"plotType"), "  Plot Type", choices = NULL,width = NULL)
            ),
            column(5,style='padding-left:3px; padding-right:12px;',align="center",
              pickerInput(NS(id,"partitionType"), "Partition", choices = NULL)
            )
          )
        ),
        uiOutput(NS(id,"box_DT")) %>% withLoader(type='html',loader = 'loader6')
      ),
      column(8,
        tabsetPanel(
          id = NS(id,"switcher"),
          type = "hidden",
          selected = "panel1",
          tabPanelBody("panel1",
            box(title = "Scatter Plot",
                width = NULL,solidHeader = T,collapsible = F,
                footer = tagList(shiny::icon("cat"), "Nya"),
              plotlyOutput(NS(id,"plot"),height = "80vh") %>% withLoader(type='html',loader = 'dnaspin')
            )
          ),
          tabPanelBody("panel2",
            box(title = "Scatter Plot",
                width = NULL,solidHeader = T,collapsible = T,
                footer = tagList(shiny::icon("cat"), "Nya"),
              dropdownButton(
                plotlyOutput(NS(id,"plot1")) %>% withLoader(type='html',loader = 'dnaspin'),
                  circle = TRUE, status = "danger", icon = icon("magnifying-glass"), width = "300px",
                  tooltip = tooltipOptions(title = "Click to choose another cluster",
                                           placement= "right"),
                  right = F
              ),
              plotlyOutput(NS(id,"plot2"),height = "80vh") %>% withLoader(type='html',loader = 'dnaspin'),
              uiOutput(NS(id,"Violin.Bar_Input")),
              conditionalPanel("typeof output.plot2 !== 'undefined'", ns = NS(id),
                               tabsetPanel(id = NS(id,"switcher2"),
                                           type = "hidden",
                                           selected = "Violin_panel",
                                           tabPanelBody("Violin_panel",
                                                        dropdownButton(
															fluidRow(
																column(7, style='padding-left:6px; padding-right:3px;',
																	column(6, style='padding-left:2px; padding-right:1px;', numericInput(NS(id,"pdf_width_violin"),"Width",value = 7)),
																	column(6, style='padding-left:1px; padding-right:2px;', numericInput(NS(id,"pdf_height_violin"),"Height",value = 7))),
																column(5, style='padding-left:0px; padding-right:6px; padding:16px', downloadButton(NS(id,'export_violin')))
															),
                                                          circle = FALSE,
                                                          status = "primary",
                                                          icon = icon("download"),
                                                          width = "300px",
                                                          size= "sm",
                                                          up = T,
                                                          tooltip = tooltipOptions(title = "Download")
                                                        ),
                                                        plotOutput(NS(id,"plot_Violin")) %>% withSpinner()
                                           ),
                                           tabPanelBody("SpikePlot_panel",
                                                        dropdownButton(
															fluidRow(
																column(7, style='padding-left:6px; padding-right:3px;',
																	column(6, style='padding-left:2px; padding-right:1px;', numericInput(NS(id,"pdf_width_SpikePlot"),"Width",value = 7)),
																	column(6, style='padding-left:1px; padding-right:2px;', numericInput(NS(id,"pdf_height_SpikePlot"),"Height",value = 7))),
																column(5, style='padding-left:0px; padding-right:6px; padding:16px', downloadButton(NS(id,'export_SpikePlot')))
															),
                                                          circle = FALSE,
                                                          status = "primary",
                                                          icon = icon("download"),
                                                          width = "300px",
                                                          size= "sm",
                                                          up = T,
                                                          tooltip = tooltipOptions(title = "Download")
                                                        ),
                                                        plotOutput(NS(id,"plot_SpikePlot")) %>% withSpinner()
                                           )
                               )
              )
              
            )
          )
        )
      )
    )
  )
}

####  Scientific notation DT ----
# js <- c(
#   "function(row, data, displayNum, index){",
#   "  var x = data[3];",
#   "  $('td:eq(3)', row).html(x.toExponential(2));",
#   "  var y = data[4];",
#   "  $('td:eq(4)', row).html(y.toExponential(2));",
#   "}"
# )

##### Marker Server Module ----
markersServer <- function(id = "markers",sce,sce.markers,point.size = 20) {
  moduleServer(id, function(input,output,session) {
    
    ### Observe Events ----
    
    updatePickerInput(session, 'partitionType', 
                      choices = names(sce.markers)
    )
    
    observeEvent(input$DTMarkers_rows_selected,{
      updateTabsetPanel(inputId = "switcher", selected = "panel2")
    })
    
    #When I change the partition if I had selected a cluster, it deleted the previous selection.
    observeEvent(c(input$partitionType,input$resetButton),ignoreInit = TRUE,{
      req(cluster_selected())
      updateTabsetPanel(inputId = "switcher", selected = "panel1")
      runjs("Shiny.setInputValue('{NS(id)('plotly_click-PlotMix')}', null);")
    })
    
    dimVector <- reactive({
      sapply(reducedDims(x = sce),FUN = ncol) 
    })
    
    observeEvent(c(dimVector()),{
      if(any(dimVector() > 3)){
        opt <- c('3','2')
      } else{
        opt <- (c("3","2")[c(3,2) %in% dimVector()])
      }
      updatePickerInput(session,inputId = "DimType", choices = opt)
    })
    
    
    observeEvent(c(input$DimType,dimVector()), {
      req(!is.null(dimVector()))
      req(input$DimType)
	  updatePickerInput(session,inputId = "plotType", choices = rev(names(which(dimVector() == as.numeric(input$DimType)  | dimVector() > 3))))
      # if(input$DimType == "3"){
        # updatePickerInput(session,inputId = "plotType", choices = rev(names(which(dimVector() == as.numeric(input$DimType)  | dimVector() > 3))))
      # } else if(input$DimType == "2") { 
        # updatePickerInput(session,inputId = "plotType", choices = rev(names(which(dimVector() == as.numeric(input$DimType)  | dimVector() > 3))))
      # }
    })
    
    observeEvent(input$Cell_Exp,{
      req(input$Cell_Exp)
      switch(input$Cell_Exp,
             'Violin'    = updateTabsetPanel(inputId = "switcher2",  selected = "Violin_panel"),
             'SpikePlot' = updateTabsetPanel(inputId = "switcher2", selected = "SpikePlot_panel")
      )
    })
    
    ### Cluster and Gen selected ----
    
    output$box_DT <- renderUI({
      if(!is.null(cluster_selected())){
        tagList(
          box(title="Cluster marker list", width = NULL, solidHeader = F,collapsible = T,
			fluidRow(column=12, align="center", style='padding-left:12px; padding-right:12px;', "Click on a gene to see its expression plots"),
			br(),
              DTOutput(NS(id,"DTMarkers"))
          ),
          fluidRow(column=12,align = "right",style='padding-left:12px; padding-right:12px;',
            actionBttn(
              inputId = NS(id,"resetButton"),
              label = "Reset", 
              style = "stretch",
              color = "primary"
            )
          )
        )
      } 
      else{
        HTML('<p style="text-align: center;"><strong>Click on a cluster to analyze its markers</strong></p>')
      }
    })
    
    cluster_selected <- reactive({
      req(!is.null(input$plotType))
      req(length(input$partitionType)>0) # same req as ClusterPLot to avoid plotly click warnings
      if(!is.null(event_data("plotly_click",source = "PlotMix"))){
        event_data("plotly_click",source = "PlotMix")[,"customdata"]
      } else{ NULL }
    })
    
    #DT markers of cluster
    output$DTMarkers <- renderDT(server = FALSE,datatable({
      req(!is.null(cluster_selected()))
      # sce.markers[[input$partitionType]][[cluster_selected()]][,c("robustness","boxcor","selectivity","meanX")]
      df <- sce.markers[[input$partitionType]][[cluster_selected()]]
      #Check if it is a cluster without any markers, if it is create a null data.frame
      if(is.null(df)){
        df <- data.frame(matrix(ncol = 3, nrow = 0))
        colnames(df) <- c('summary.stats','log.FDR','boxcor')
      }
      df
      
    }, selection = 'single',
    extensions = 'Buttons',
    options = list(language = list(zeroRecords = "No markers found for this cluster or criterion"),
                   dom = 'Bfrtip',
                   exportOptions = list(header = ""),
                   buttons = c('copy', 'csv', 'excel', 'pdf'),
                   # rowCallback = DT::JS(js),
                   scrollX=T,
                   autoWidth = F),
    rownames= TRUE,
    caption = htmltools::tags$caption(
      style = paste0('caption-side: top; text-align: center; font-weight: bold;color:white;background-color:',OrderPartReact()$colPart[[cluster_selected()]]),
      cluster_selected())
    ) %>% formatRound(columns = c('boxcor'),digits = 3)
    )
    
    #Gene selected from the DT
    gene_marker_selected <- eventReactive(input$DTMarkers_rows_selected,{
      req(!is.null(input$DTMarkers_rows_selected))
      rownames(sce.markers[[input$partitionType]][[cluster_selected()]][input$DTMarkers_rows_selected,])
    })
    
    ###Plots ----
    
    OrderPartReact <- eventReactive(input$partitionType,{
      req(input$partitionType)
      Col.and.Order(partition = input$partitionType, sce=sce)
    })
    
          #### Scatter ----
    ClusterPlot <- eventReactive(c(input$DimType,input$plotType,input$partitionType),{
      req(!is.null(input$plotType))
      req(length(input$partitionType)>0)
      #3D
      if(input$DimType == "3"){
        plot_ly(type = "scatter3d", mode = "markers",source = "PlotMix")  %>%
          layout(dragmode = "select",
                 scene = list(xaxis = list(title = 'Dim1',showgrid=F,visible=F),
                               yaxis = list(title = 'Dim2',showgrid=F,visible=F),
                               zaxis = list(title = 'Dim3',showgrid=F,visible=F)),
                 title = paste(input$partitionType, 'Partition'),
                 margin = list(l = 0,
                               r = 10,
                               b = 0,
                               t = 40,
                               pad = 0)) %>% 
          add_markers(x = ~reducedDim(sce,input$plotType)[,1], y=~reducedDim(sce,input$plotType)[,2], z=~reducedDim(sce,input$plotType)[,3],
                      color = ~I(as.character(OrderPartReact()$colPart[colData(sce)[,input$partitionType]])),
                      name = ~colData(sce)[,input$partitionType],
                      customdata= ~colData(sce)[,input$partitionType],
                      size = I(point.size),span=I(0),text=~colData(sce)[,input$partitionType],hoverinfo='text') %>% 
          event_register("plotly_selecting") %>% 
          toWebGL()
      } else { #2D
        req(!is.null(input$plotType))
        req(!is.null(input$partitionType))
        plot_ly(type = "scatter", mode = "markers",source="PlotMix")  %>%
          layout(dragmode = "select",
                 xaxis = list(title = 'Dim1',zeroline=F),
                 yaxis = list(title = 'Dim2',zeroline=F),
                 title = paste(input$partitionType, 'Partition')) %>%
          add_markers(x = ~reducedDim(sce,input$plotType)[,1], y=~reducedDim(sce,input$plotType)[,2], 
                      color = ~I(as.character(OrderPartReact()$colPart[colData(sce)[,input$partitionType]])),
                      name = ~colData(sce)[,input$partitionType],
                      customdata= ~colData(sce)[,input$partitionType], 
                      size = I(point.size),span=I(0),text=~colData(sce)[,input$partitionType],hoverinfo='text') %>% 
          event_register("plotly_selecting") %>%
		  config(modeBarButtonsToRemove = c("select2d", "lasso2d", "hoverCompareCartesian")) %>%
          toWebGL()
      }
      
    })
    
    ExpressionPlot <- eventReactive(c(input$DimType,input$plotType,input$partitionType,input$DTMarkers_rows_selected),{
      req(input$DTMarkers_rows_selected)
      #3D
      if(input$DimType == "3"){
        plot_ly(type = "scatter3d", mode = "markers")  %>%
          layout(dragmode = "select",
                 scene = list(xaxis = list(title = 'Dim1',showgrid=F,visible=F),
                               yaxis = list(title = 'Dim2',showgrid=F,visible=F),
                               zaxis = list(title = 'Dim3',showgrid=F,visible=F)),
                 legend= list(x=1,y=1),
                 showlegend = TRUE,
                 title = paste(gene_marker_selected(), 'Expression'),
                 margin = list(l = 0,
                               r = 10,
                               b = 0,
                               t = 40,
                               pad = 0)) %>% 
          add_markers(x = ~reducedDim(sce,input$plotType)[,1], y=~reducedDim(sce,input$plotType)[,2], 
                      z=~reducedDim(sce,input$plotType)[,3],
                      text=~colData(sce)[,input$partitionType],
                      hoverinfo='text',
                      customdata= ~colData(sce)[,input$partitionType],
                      color = ~logcounts(sce)[gene_marker_selected(),],
                      size = I(point.size),span=I(0),
                      name = ~colData(sce)[,input$partitionType]) %>% 
          colorbar(title = "log(counts)",x=0,y=1) %>% 
          toWebGL()
      } else { #2D
        plot_ly(type = "scatter", mode = "markers")  %>%
          layout(dragmode = "select",
                 xaxis = list(title = 'Dim1',zeroline=F),
                 yaxis = list(title = 'Dim2',zeroline=F),
                 showlegend = FALSE,
                 title = paste(gene_marker_selected(), 'Expression')) %>% 
          add_markers(x = ~reducedDim(sce,input$plotType)[,1], y=~reducedDim(sce,input$plotType)[,2],
                      customdata= ~colData(sce)[,input$partitionType],
                      color = ~logcounts(sce)[gene_marker_selected(),], 
                      size = I(point.size),span=I(0),text=~colData(sce)[,input$partitionType],hoverinfo='text') %>% 
          colorbar(title = "log(counts)") %>%
		  config(modeBarButtonsToRemove = c("select2d", "lasso2d", "hoverCompareCartesian")) %>% 
          toWebGL()
      }
      
    })
    
    output$plot <- renderPlotly(ClusterPlot())
    
    output$plot1 <- renderPlotly({
      ClusterPlot() %>% layout(showlegend = FALSE)})
    
    output$plot2 <- renderPlotly({
      req(input$DTMarkers_rows_selected)
      ExpressionPlot()
    })
    
    
          ####  Violin&SpikePlots ----
    
    output$Violin.Bar_Input <- renderUI({
      req(!is.null(gene_marker_selected()))
      radioGroupButtons(inputId = NS(id,"Cell_Exp"), label=NULL,choices = c("Violin", "SpikePlot"),
                        direction = "horizontal",justified = T,individual=T)
    })
    
    ViolinReact <- reactive({
      req(!is.null(gene_marker_selected()))
      data.frame(Y = logcounts(sce)[gene_marker_selected(),],
                 X=factor(colData(sce)[,input$partitionType]))
    })
    
    ViolinReact_Cell <- reactive({
      ViolinReact() %>% group_by(X) %>% summarise(n=n(),Ymax = (max(Y)+0.5))
    })
    
    ViolinPlot <-reactive({
      req(!is.null(gene_marker_selected()))
      ggplot(ViolinReact()) + 
        geom_violin(aes(y = Y, 
                        x = X, 
                        fill = X),
                    scale="width") + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              legend.position = "none") + 
        xlab("Clusters") + 
        ylab("log(counts)") +
        ggtitle(paste(gene_marker_selected(), 'expression')) +
        geom_text(aes(label = n,x=X, y=Ymax), data = ViolinReact_Cell()) +
        scale_fill_manual(values=OrderPartReact()$colPart)
    })
    
    SpikePlot <-reactive({
      req(!is.null(gene_marker_selected()))
      m <- barplot(assay(sce,"logcounts")[gene_marker_selected(),][OrderPartReact()$ordPart],
              col = OrderPartReact()$colPart[colData(sce)[,input$partitionType]][OrderPartReact()$ordPart],
              border = OrderPartReact()$colPart[colData(sce)[,input$partitionType]][OrderPartReact()$ordPart],
              ylab = "log(counts)", main = paste(gene_marker_selected(), 'expression'), names.arg = F)
	  colleg <- legend_col(names(OrderPartReact()$colPart), max(m))
      legend(max(m)/2, -0.05, legend = names(OrderPartReact()$colPart), col = OrderPartReact()$colPart,
             pch=19, xpd=T, xjust = 0.5, cex = 0.9, ncol=colleg$ncol, text.width = colleg$colwidth)
      lines(x = m,
            tapply(assay(sce,"logcounts")[gene_marker_selected(),],
                   INDEX = colData(sce)[,input$partitionType],
                   FUN = mean)[colData(sce)[,input$partitionType]][OrderPartReact()$ordPart],
            lty=2,col="black")
      graph <- recordPlot()
      graph
    })
    
    output$plot_Violin <- renderPlot({
      req(!is.null(ViolinPlot()))
      ViolinPlot() %>% plot()
    })
    
    output$plot_SpikePlot <- renderPlot({
      req(!is.null(SpikePlot()))
      SpikePlot() %>% print()
    })
    
    ### Downloads -----
    
    output$export_violin = downloadHandler(
      filename = function() {"Violin_Markers.pdf"},
      content = function(file) {
        pdf(file,
            width = input$pdf_width_violin,
            height = input$pdf_height_violin
        )
        ViolinPlot() %>% plot()
        dev.off()
      })
    
    output$export_SpikePlot = downloadHandler(
      filename = function() {"SpikePlot_Markers.pdf"},
      content = function(file) {
        pdf(file,
            width = input$pdf_width_SpikePlot,
            height = input$pdf_height_SpikePlot
        )
        SpikePlot() %>% print()
        dev.off()
      })
  })
}

