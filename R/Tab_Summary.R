#####                              ######## 
#####     QualityControl Tab       ######## 
#####                              ######## 

##### QC UI Module ----
QC_UI <- function(id) {
  tagList(
    uiOutput(outputId = NS(id,"DescriptionText")),
    box(width = 12,title = "Summary",solidHeader = F, status = "primary",collapsible = T,
      fluidRow(
        infoBoxOutput(NS(id,"CellBox")) %>% withSpinner(),
        infoBoxOutput(NS(id,"GenesBox")),
        infoBoxOutput(NS(id,"nFeatureBox")),
      ),
      fluidRow(
        infoBoxOutput(NS(id,"PartitionBox")) ,
        infoBoxOutput(NS(id,"ExpBox")),
        infoBoxOutput(NS(id,"nGeneBox"))
      )
    ),
    fluidRow(
      column(2,
        box(title = htmltools::span(icon("gears"), " Settings"),
            width = NULL, status = "primary",solidHeader = T,collapsible = T,
          pickerInput(
            inputId = NS(id,"partitionType"),
            label = "Partition", 
            choices = NULL
          ),
          radioGroupButtons(inputId = NS(id,"scatter_box"),
                            label = "",
                            choices = c("Scatter", "Boxplot"),
                            status = "primary",
                            direction = "horizontal",justified = T)
        )
      ),
      column(10,
        box( width = NULL, title = "Quality Control",
             solidHeader = T,collapsible = T,
             footer = tagList(shiny::icon("cat"), "Nya"),
             tabsetPanel(id = NS(id,"switcher2"),
                         type = "hidden",
                         selected = "Scatter_panel",
                         tabPanelBody("Scatter_panel",
                                      dropdownButton(
										fluidRow(
											column(7, style='padding-left:6px; padding-right:3px;',
												column(6, style='padding-left:2px; padding-right:1px;', numericInput(NS(id,"pdf_width_scatter"),"Width",value = 7)),
												column(6, style='padding-left:1px; padding-right:2px;', numericInput(NS(id,"pdf_height_scatter"),"Height",value = 7))),
											column(5, style='padding-left:0px; padding-right:6px; padding:16px;',downloadButton(NS(id,'export_scatter')))
										),
                                        circle = FALSE,
                                        status = "primary",
                                        icon = icon("download"),
                                        width = "300px",
                                        size= "sm",
                                        up = F,
                                        tooltip = tooltipOptions(title = "Download")
                                      ),
                                      plotlyOutput(NS(id,"plot_scatter"),
                                                   height = "80vh") %>% withLoader(type='html',loader = 'dnaspin')
                         ),
                         tabPanelBody("boxPlot_panel",
                                      dropdownButton(
										fluidRow(
											column(7, style='padding-left:6px; padding-right:3px;',
												column(6, style='padding-left:2px; padding-right:1px;', numericInput(NS(id,"pdf_width_boxplot"),"Width",value = 7)),
												column(6, style='padding-left:1px; padding-right:2px;', numericInput(NS(id,"pdf_height_boxplot"),"Height",value = 7))),
											column(5, style='padding-left:0px; padding-right:6px; padding:16px;', downloadButton(NS(id,'export_boxplot')))
										),
                                        circle = FALSE,
                                        status = "primary",
                                        icon = icon("download"),
                                        width = "300px",
                                        size= "sm",
                                        up = F,
                                        tooltip = tooltipOptions(title = "Download")
                                      ),
                                      plotlyOutput(NS(id,"plot_boxPlot"),
                                                   height = "80vh") %>% withLoader(type='html',loader = 'dnaspin')
                         )
             )
        )
      )
    )
  )
}

##### Expression Server Module ----
QC_Server <- function(id,sce,descriptionText) {
  moduleServer(id, function(input,output,session) {
    
    ### Observe Events ----
    updatePickerInput(session, 'partitionType', 
                      choices = c("None",names(colData(sce))[sapply(colData(sce), is.factor)]))
    
    observeEvent(input$scatter_box,{
      req(input$scatter_box)
      switch(input$scatter_box,
             'Scatter'    = updateTabsetPanel(inputId = "switcher2",  selected = "Scatter_panel"),
             'Boxplot'    = updateTabsetPanel(inputId = "switcher2", selected = "boxPlot_panel")
      )
    })
    
    ### Description Text -----
    output$DescriptionText <- renderUI({
      req(!is.null(descriptionText))
      box(width = 12,title = "Description",solidHeader = T, status = "primary",collapsible = T,
          p(descriptionText)
      )
    })
    
    #### Summary Box ----
    output$CellBox <- renderInfoBox({
      infoBox("Cells", value = ncol(sce),
              icon = icon('vial'), #icon("sign-in"),
              color = "light-blue", fill = F
      )
    })
    
    output$GenesBox <- renderInfoBox({
      infoBox(title = "Genes", value = nrow(sce),
              icon = icon('dna'), #icon("sign-in"),
              color = "light-blue", fill = F
      )
    })
    
    output$PartitionBox <- renderInfoBox({
      infoBox(title = "Partitions", value = ncol(colData(sce)[sapply(colData(sce), is.factor)]),
              icon = icon('filter'), #icon("sign-in"),
              color = "light-blue", fill = F
      )
    })
    
    output$ExpBox <- renderInfoBox({
      infoBox(title="Experiments",value= length(sce@metadata),subtitle = "Integrated experiments",
              icon = icon('coins'), #icon("sign-in"),
              color = "light-blue", fill = F
      )
    })
    
    output$nFeatureBox <- renderInfoBox({
      infoBox(title="Gene Detected",value= round(mean(sce$nFeatures),2),subtitle = "Mean gene detected per cell",
              icon = icon('magnifying-glass-chart'),
              color = "light-blue", fill = F
      )
    })
    
    output$nGeneBox <- renderInfoBox({
      infoBox(title="Lib size",value= round(mean(sce$nCounts),2),subtitle = "Mean lib size per cell",
              icon = icon('book'), #icon("sign-in"),
              color = "light-blue", fill = F
      )
    })
    
    #### Plots ----
    OrderPartReact <- eventReactive(input$partitionType,{
      if(input$partitionType == "None"){
        list(colPart = "grey")
      } else{
        Col.and.Order(partition = input$partitionType, sce=sce)
      }
    })
    
    #### Scatter & Boxplot ----
    ScatterPlot <- reactive({
        by_color <- if(input$partitionType == "None") {NULL} else{colData(sce)[,input$partitionType]}
        g <- ggplot(as.data.frame(colData(sce)),aes(x=nCounts,y=nFeatures)) + geom_point(alpha=0.5,aes(col=by_color))+ 
             xlab('nCounts') + ylab('nFeatures') + 
             labs(color = input$partitionType) + 
             scale_colour_manual(values=OrderPartReact()$colPart)
          
        g
      })
    
    output$plot_scatter <- renderPlotly({
      req(!is.null(ScatterPlot()))
      g <- ggplotly(ScatterPlot()) %>% config(modeBarButtonsToRemove = c("select2d", "lasso2d"))
      g %>% toWebGL()
    })
    
    BoxPlot <- reactive({
      
        df <-data.frame(Yratio=sce$nCounts/sce$nFeatures, Counts = sce$nCounts,Features  = sce$nFeatures)
        if(input$partitionType == "None") {
          df$X <- "All"
        } else {
          df$X <- colData(sce)[,input$partitionType]
        }
        
        summ <- df  %>% group_by(X) %>% summarise(n=n(),Ymax = (max(Yratio)+1))
        
        g_ratio <- ggplot(df, aes(y=Yratio,x=X,fill=X)) + geom_boxplot() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none") + geom_text(aes(label = n,y=Ymax), data = summ) +
          ylab('nCounts/nFeatures') +
          xlab(input$partitionType) + scale_fill_manual(values=OrderPartReact()$colPart)
        
        # g_ratio <- ggplotly(g_ratio)
        
        g_counts <- ggplot(df, aes(y=Counts,x=X,fill=X)) + geom_boxplot() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none") +
          ylab('nCounts') +
          xlab(input$partitionType) + scale_fill_manual(values=OrderPartReact()$colPart)
        # g_counts <- g_counts %>% ggplotly()
        
        g_feat <- ggplot(df, aes(y=Features,x=X,fill=X)) + geom_boxplot() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none") + 
          ylab('nFeatures') +
          xlab(input$partitionType) + scale_fill_manual(values=OrderPartReact()$colPart)
        
        # g_feat <- g_feat %>% ggplotly()
        
        
      list(g1 = g_ratio,g2= g_counts,g3= g_feat)        
       #  g <- subplot(g_ratio,g_counts,g_feat,nrows = 3,shareX = T,shareY = F,titleY = T,heights = c(0.4,0.3,0.3)) %>% 
       #    config(modeBarButtonsToRemove = c("select2d", "lasso2d"))
       #  
       # g %>% toWebGL()
       
    })
    
    output$plot_boxPlot <- renderPlotly({
      req(!is.null(BoxPlot()))
      g <- subplot(ggplotly(BoxPlot()$g1),ggplotly(BoxPlot()$g2),ggplotly(BoxPlot()$g3),nrows = 3,shareX = T,shareY = F,titleY = T,heights = c(0.4,0.3,0.3)) %>% 
        config(modeBarButtonsToRemove = c("select2d", "lasso2d"))
      g %>% toWebGL()
    })
    
    #### Download ----
    output$export_scatter = downloadHandler(
      filename = function() {"ScatterPlot_QC.pdf"},
      content = function(file) {
        pdf(file,
            width = input$pdf_width_scatter,
            height = input$pdf_height_scatter
        )
        ScatterPlot() %>% plot()
        dev.off()
      }
    )
    
    output$export_boxplot = downloadHandler(
      filename = function() {"BoxPlot_QC.pdf"},
      content = function(file) {
        pdf(file,
            width = input$pdf_width_boxplot,
            height = input$pdf_height_boxplot
        )
        
        g1 <- BoxPlot()$g1 + theme(
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()
        )

        g2 <- BoxPlot()$g2 + theme(
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()
        )
        
        g <- cowplot::plot_grid(
          g1, g2, BoxPlot()$g3,
          ncol = 1,  # Number of columns in the grid
          align = "v",  # Align the plots vertically
          axis = "x",  # Display left and right axes
          nrow = 3,  # Number of rows in the grid
          rel_heights = c(0.3,0.3,0.4),  # Heights for each row
          labels = NULL  # Automatically label the plots
        )
        
        g %>% plot()
        dev.off()
      }
    )
  })
}
