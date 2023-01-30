#####                              ######## 
#####        Markers Tab           ######## 
#####                              ######## 

##### Marker UI Module ----
markersUI <- function(id = "markers") {
  tagList(
    useShinyjs(),
    fluidRow(
      column(4,
        box(title = htmltools::span(icon("fa-light fa-gears"), " Settings"),
            width = NULL, status = "primary",solidHeader = T,collapsible = T,
          fluidRow(
            column(4,style='padding-left:12px; padding-right:3px;',align="center",
              pickerInput(NS(id,"plotType"), "  Plot Type", choices = NULL,width = NULL)
            ),
            column(3,style='padding-left:3px; padding-right:3px;',align="center",
              pickerInput(NS(id,"DimType"), "  Dim Type", choices = NULL,width = NULL)
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
            box(title = "Scatter",
                width = NULL,solidHeader = T,collapsible = F,
                footer = tagList(shiny::icon("cat"), "Nya"),
              plotlyOutput(NS(id,"plot"),height = "100vh") %>% withLoader(type='html',loader = 'dnaspin')
            )
          ),
          tabPanelBody("panel2",
            box(title = "Scatter",
                width = NULL,solidHeader = T,collapsible = T,
                footer = tagList(shiny::icon("cat"), "Nya"),
              dropdownButton(
                plotlyOutput(NS(id,"plot1")) %>% withLoader(type='html',loader = 'dnaspin'),
                  circle = TRUE, status = "danger", icon = icon("fa-solid fa-magnifying-glass"), width = "300px",
                  tooltip = tooltipOptions(title = "Click to choose another cluster",
                                           placement= "right"),
                  right = F
              ),
              plotlyOutput(NS(id,"plot2"),height = "100vh") %>% withLoader(type='html',loader = 'dnaspin')
            ),
            box(title = "Expression",
                width = NULL,solidHeader = T,collapsible = T,
                footer = tagList(shiny::icon("cat"), "Nya"),
              uiOutput(NS(id,"Violin.Bar_Input")),
              plotOutput(NS(id,"Violin.Bar_Plot")) %>%  shinycssloaders::withSpinner()
            )
          )
        )
      )
    )
  )
}

##### Marker Server Module ----
markersServer <- function(id = "markers",sce,ldf,point.size = 20) {
  moduleServer(id, function(input,output,session) {
    
    ### Observe Events ----
    
    updatePickerInput(session, 'partitionType', 
                      choices = names(colData(sce))[sapply(colData(sce), is.factor)])
    
    observeEvent(input$DTMarkers_rows_selected,{
      updateTabsetPanel(inputId = "switcher", selected = "panel2")
    })
    
    #When I change the partition if I had selected a cluster, it deleted the previous selection.
    observeEvent(c(input$partitionType,input$resetButton),ignoreInit = TRUE,{
      req(cluster_selected())
      updateTabsetPanel(inputId = "switcher", selected = "panel1")
      runjs("Shiny.setInputValue('plotly_click-PlotMix', null);")
    })
    
    dimVector <- reactive({
      sapply(reducedDims(x = sce),FUN = ncol) 
    })
    
    observeEvent(c(dimVector()),{
      updatePickerInput(session,inputId = "DimType", choices = (c("3","2")[c(3,2) %in% dimVector()]))
    })
    
    
    observeEvent(c(input$DimType,dimVector()), {
      req(!is.null(dimVector()))
      req(input$DimType)
      if(input$DimType == "3"){
        updatePickerInput(session,inputId = "plotType", choices = rev(names(which(dimVector() == as.numeric(input$DimType)  | dimVector() > 3))))
      } else if(input$DimType == "2") { 
        updatePickerInput(session,inputId = "plotType", choices = rev(names(which(dimVector() == as.numeric(input$DimType)  | dimVector() > 3))))
      }
    })
    
    ### Cluster and Gen selected ----
    
    output$box_DT <- renderUI({
      if(!is.null(cluster_selected())){
      tagList(
      box(width = NULL, status = "primary",solidHeader = F,collapsible = T,
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
      )} else{
        
        HTML('<p style="text-align: center;"><strong>Click on a cluster to analyze its markers</strong></p>')
      }
    })
    
    cluster_selected <- reactive({
      req(input$partitionType)
      if(!is.null(event_data("plotly_click",source = "PlotMix"))){
        event_data("plotly_click",source = "PlotMix")[,"customdata"]
      } else{ NULL }
    })
    
    #DT markers of cluster
    output$DTMarkers <- renderDT(server = FALSE,datatable({
      req(!is.null(cluster_selected()))
      # ldf[[input$partitionType]][[cluster_selected()]][,c("robustness","boxcor","selectivity","meanX")]
      ldf[[input$partitionType]][[cluster_selected()]]
    }, selection = 'single',
    extensions = 'Buttons',
    options = list(language = list(zeroRecords = "Click on a cluster to analyze its markers (double-click to clear)"),
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
    ) %>% formatRound(columns = c("boxcor"),digits = 3)
    )
    
    #Gene selected from the DT
    gene_marker_selected <- eventReactive(input$DTMarkers_rows_selected,{
      req(!is.null(input$DTMarkers_rows_selected))
      rownames(ldf[[input$partitionType]][[cluster_selected()]][input$DTMarkers_rows_selected,])
    })
    
    ###Plots ----
    
    OrderPartReact <- eventReactive(input$partitionType,{
      req(input$partitionType)
      Col.and.Order(partition = input$partitionType, sce=sce)
    })
    
          #### Scatter ----
    ClusterPlot <- eventReactive(c(input$plotType,input$partitionType),{
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
          toWebGL()
      }
      
    })
    
    ExpressionPlot <- eventReactive(c(input$plotType,input$partitionType,input$DTMarkers_rows_selected),{
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
      req(!is.null(gene_marker_selected()) & input$Cell_Exp == "Violin")
      data.frame(Y = logcounts(sce)[gene_marker_selected(),],
                 X=factor(colData(sce)[,input$partitionType]))
    })
    
    ViolinReact_Cell <- reactive({
      ViolinReact() %>% group_by(X) %>% summarise(n=n(),Ymax = (max(Y)+0.5))
    })
    
    ViolinPlot <-reactive({
      req(!is.null(gene_marker_selected()) & input$Cell_Exp == "Violin")
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
      req(!is.null(gene_marker_selected()) & input$Cell_Exp == "SpikePlot")
      barplot(assay(sce,"logcounts")[gene_marker_selected(),][OrderPartReact()$ordPart],
              col = OrderPartReact()$colPart[colData(sce)[,input$partitionType]][OrderPartReact()$ordPart],
              border = OrderPartReact()$colPart[colData(sce)[,input$partitionType]][OrderPartReact()$ordPart],
              ylab = "log(counts)", main = paste(gene_marker_selected(), 'expression'), names.arg = F) 
      legend("bottom", legend = names(OrderPartReact()$colPart), col = OrderPartReact()$colPart,
             pch=19, ncol=6, xpd=T, inset=c(0,-0.25))
      abline(h=mean(assay(sce,"logcounts")[gene_marker_selected(),][assay(sce,"logcounts")[gene_marker_selected(),]>0]),lty=2,col="grey")
    })
    
    output$Violin.Bar_Plot <- renderPlot({
      req(!is.null(gene_marker_selected()))
      req(input$Cell_Exp)
      if(input$Cell_Exp == "Violin"){
        ViolinPlot()
      } else { 
        SpikePlot()
      }
    })
    
    
  })
}

