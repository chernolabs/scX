#####                              ######## 
#####       NewMarkers Tab         ######## 
#####                              ######## 

##### New Marker UI Module ----
N_markersUI <- function(id) {
  tagList(
    useShinyjs(),
    fluidRow(
      column(3,
        box(title = htmltools::span(icon("gears"), " Settings"), 
            width = NULL, status = "primary",solidHeader = T,collapsible = T,
          fluidRow(
            column(6,style='padding-left:12px; padding-right:3px;', align="center",
              pickerInput(NS(id,"plotType2D"), 
                          "Plot Type 2D", 
                          choices = NULL)
            ),
            column(6,style='padding-left:3px; padding-right:12px;', align="center",
              pickerInput(NS(id,"partitionType"),
                          "Partition",
                          choices = NULL)
            )
          ),
          conditionalPanel(ns = NS(id), "typeof output.plot2 !== 'undefined' && input.switcher == 'panel2'",
			hr(style = "border-top: 1px solid #0073b7;"),
            fluidRow(style='padding-left:12px; padding-right:12px;', h4("Gene expression scatter plot")),
			fluidRow(
              column(6,style='padding-left:12px; padding-right:3px;', align="center",
                pickerInput(NS(id,"DimType"),
                            "  # dims",
                            choices = NULL,
                            width = NULL)
              ),
              column(6,style='padding-left:3px; padding-right:12px;', align="center",
                pickerInput(NS(id,"plotType"), 
                            "  Plot Type", 
                            choices = NULL,
                            width = NULL)
              )
            )
          ),
          conditionalPanel(ns = NS(id), "typeof output.DTMarkers !== 'undefined'", 
                           fluidRow(
                             column(12,align ="center",style='padding-left:12px; padding-right:12px;',
                                    downloadButton(NS(id,'export'),label = "Download Selected Cells")
                             )
                           )
          )
        ), 
        uiOutput(NS(id,"box_DT")) %>% withLoader(type='html',loader = 'loader6')
      ),
      column(9,
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
                tooltip = tooltipOptions(title = "Click to select other group of cells",
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
																column(4, style='padding-left:12px; padding-right:3px;', numericInput(NS(id,"pdf_width_violin"),"Width",value = 7)),
																column(4, style='padding-left:3px; padding-right:3px;', numericInput(NS(id,"pdf_height_violin"),"Height",value = 7)),
																column(4, style='padding-left:3px; padding-right:12px; padding:16px', downloadButton(NS(id,'export_violin')))
															),
                                                          circle = FALSE,
                                                          status = "primary",
                                                          icon = icon("download"),
                                                          width = "300px",
                                                          size= "sm",
                                                          up = T,
                                                          tooltip = tooltipOptions(title = "Press to Download")
                                                        ),
                                                        plotOutput(NS(id,"plot_Violin")) %>% withSpinner()
                                           ),
                                           tabPanelBody("SpikePlot_panel",
                                                        dropdownButton(
															fluidRow(
																column(4, style='padding-left:12px; padding-right:3px;', numericInput(NS(id,"pdf_width_SpikePlot"),"Width",value = 7)),
																column(4, style='padding-left:3px; padding-right:3px;', numericInput(NS(id,"pdf_height_SpikePlot"),"Height",value = 7)),
																column(4, style='padding-left:3px; padding-right:12px; padding:16px', downloadButton(NS(id,'export_SpikePlot')))
															),
                                                          circle = FALSE,
                                                          status = "primary",
                                                          icon = icon("download"),
                                                          width = "300px",
                                                          size= "sm",
                                                          up = T,
                                                          tooltip = tooltipOptions(title = "Press to Download")
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

##### New Marker UI Module ----
N_markersServer <- function(id,sce,point.size = 20) {
  moduleServer(id, function(input,output,session) {
    
    ### Observe Events ----
    
    updatePickerInput(session, 'partitionType', 
                      choices = names(colData(sce))[sapply(colData(sce), is.factor)])
    
    
    dimVector <- reactive({
      sapply(reducedDims(x = sce),FUN = ncol) 
    })
    
    observeEvent(dimVector(),{
      updatePickerInput(session,inputId = "plotType2D", choices = rev(names(which(dimVector() == 2  | dimVector() > 3))))
    })
    
    observeEvent(c(dimVector(),input$switcher),{
      updatePickerInput(session,inputId = "DimType", choices = (c("3","2")[c(3,2) %in% dimVector()]))
    })
    
    observeEvent(c(input$DimType,input$switcher), {
      req(!is.null(dimVector()))
      req(input$DimType)
      if(input$DimType == "3"){
        updatePickerInput(session,inputId = "plotType", choices = rev(names(which(dimVector() == as.numeric(input$DimType)  | dimVector() > 3))))
      } else if(input$DimType == "2") { 
        updatePickerInput(session,inputId = "plotType", choices = rev(names(which(dimVector() == as.numeric(input$DimType)  | dimVector() > 3))))
      }
    })
    
    observeEvent(gene_marker_selected(),{
      updateTabsetPanel(inputId = "switcher", selected = "panel2")
    })
    
    #When I change the partition if I had selected a cluster, it deleted the previous selection.
    observeEvent(input$resetButton,ignoreInit = TRUE,{
      req(!is.null(cells_selected()))
      updateTabsetPanel(inputId = "switcher", selected = "panel1")
      runjs("Shiny.setInputValue('plotly_selected-PlotMix', null);")
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
      if(!is.null(MarkersDT())){
        tagList(
          box(title = "Selection marker list", width = NULL,solidHeader = F,collapsible = T,
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
      } else{
        HTML('<p style="text-align: center;"><strong>Select custom clusters using the box and lasso select tools on the scatter plot.</strong></p>')
      }
    })
    
    cells_selected <- reactive({
      req(input$partitionType)
      d <- event_data(source = "PlotMix","plotly_selected")
      # d <- d[,"customdata"]
      # if (is.null(d)) NULL else d
      if(!is.null(d)  & length(d)>0){
        d[,"customdata"]
      } else{ NULL }
    })
    
    MarkersDT <- reactive({
      if(!is.null(cells_selected())){
      a <- cajitasdeluz(ssce = sce,selected.cells =cells_selected(),corr = 0.3)  
      a
      } 
      else { NULL }
    })
    
    # #DT markers of cluster
    output$DTMarkers <- renderDT(server = FALSE,
                                 datatable({
                                   req(!is.null(MarkersDT()))
                                   MarkersDT()
                                 }, selection = 'single',
                                 extensions = 'Buttons',
                                 
                                 options = list(language = list(zeroRecords = "Click on a cluster to analyze its markers (double-click to clear)"),
                                                dom = 'Bfrtip',
                                                exportOptions = list(header = ""),
                                                buttons = c('copy', 'csv', 'excel', 'pdf')
                                                # rowCallback = JS(js)),
                                 ),
                                 rownames= TRUE
                                 ) %>% formatRound(columns = c("box.cor"),digits = 3)
    )
    # #Gene selected from the DT
    gene_marker_selected <- eventReactive(input$DTMarkers_rows_selected,{
      req(!is.null(input$DTMarkers_rows_selected))
      rownames(MarkersDT()[input$DTMarkers_rows_selected,,drop=F])
    })
     
    ### Download selected cells --
    
    output$export = downloadHandler(
      filename = function() {"Selected_Cells.csv"},
      content = function(file) {
        req(!is.null(MarkersDT()))
        write.csv(data.frame(Selected_Cells = cells_selected()),
                  file = file,
                  row.names = F)
      }
    )
    
    ###Plots ----
    
    OrderPartReact <- eventReactive(input$partitionType,{
      req(input$partitionType)
      Col.and.Order(partition = input$partitionType, sce=sce)
    })
    
          #### Scatter ----
    ClusterPlot <- eventReactive(c(input$plotType2D,input$partitionType),{
      req(!is.null(OrderPartReact()))

        plot_ly(type = "scatter", mode = "markers",source="PlotMix")  %>%
          layout(dragmode = "select",
                 xaxis = list(title = 'Dim1',zeroline=F),
                 yaxis = list(title = 'Dim2',zeroline=F),
                 title = paste(input$partitionType, 'Partition')) %>%
          add_markers(x = ~reducedDim(sce,input$plotType2D)[,1], y=~reducedDim(sce,input$plotType2D)[,2], 
                      color = ~I(as.character(OrderPartReact()$colPart[colData(sce)[,input$partitionType]])),
                      name = ~colData(sce)[,input$partitionType],
                      customdata= ~I(colnames(sce)), 
                      size = I(point.size),span=I(0),text=~colData(sce)[,input$partitionType],hoverinfo='text') %>% 
          event_register("plotly_selecting") %>% 
          toWebGL()
      
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
          add_markers(x = ~reducedDim(sce,input$plotType)[,1], y=~reducedDim(sce,input$plotType)[,2], z=~reducedDim(sce,input$plotType)[,3],
                      # customdata= ~colData(sce)[,input$partitionType],
                      name = ~colData(sce)[,input$partitionType],
                      color = ~logcounts(sce)[gene_marker_selected(),],
                      size = I(point.size),span=I(0),text=~colData(sce)[,input$partitionType],hoverinfo='text') %>% 
          colorbar(title = "log(counts)",x=0,y=1) %>% 
          toWebGL()
      } else { #2D
        plot_ly(type = "scatter", mode = "markers")  %>%
          layout(dragmode = "select",
                 xaxis = list(title = 'Dim1',zeroline=F),
                 yaxis = list(title = 'Dim2',zeroline=F),
                 legend= list(x=1,y=1),
                 showlegend = FALSE,
                 title = paste(gene_marker_selected(), 'Expression')) %>% 
          add_markers(x = ~reducedDim(sce,input$plotType)[,1], y=~reducedDim(sce,input$plotType)[,2],
                      # customdata= ~colData(sce)[,input$partitionType],
                      #name = ~colData(sce)[,input$partitionType],
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
    
          #### Violin&SpikePlots ------
    
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
      legend("bottom", legend = names(OrderPartReact()$colPart), col = OrderPartReact()$colPart,
             pch=19, ncol=6, xpd=T, inset=c(0,-0.25))
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
      filename = function() {"Violin_NewMarkers.pdf"},
      content = function(file) {
        pdf(file,
            width = input$pdf_width_violin,
            height = input$pdf_height_violin
        )
        ViolinPlot() %>% plot()
        dev.off()
      })
    
    output$export_SpikePlot = downloadHandler(
      filename = function() {"SpikePlot_NewMarkers.pdf"},
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

