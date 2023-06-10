#####                              ######## 
#####       CoExpression Tab       ######## 
#####                              ######## 

##### CoExpression UI Module ----
COExpUI <- function(id) {
  tagList(
    fluidRow(
      column(3,
        box(title = htmltools::span(icon("gears"), " Settings"),
            width = NULL, status = "primary",solidHeader = T,collapsible = F,
          fluidRow(
            column(6,style='padding-left:12px; padding-right:3px;', align="center",
              pickerInput(NS(id,"DimType"), "  # dims", choices = NULL,width = NULL)
            ),
            column(6,style='padding-left:3px; padding-right:12px;', align="center",
              pickerInput(NS(id,"plotType"), "  Plot Type", choices = NULL,width = NULL)
            )
          ),
          conditionalPanel("typeof output.plot_expression !== 'undefined'", ns = NS(id), 
            fluidRow(
              column(8,style='padding-left:12px; padding-right:3px;',
                pickerInput(inputId = NS(id,"partitionType"), 
                            label = "Partition",
                            choices = NULL)
              ),
              column(4,style='padding-left:3px; padding-right:1px;padding-top:12px',
                br(),
                prettyCheckbox(NS(id,"button"),
                               label="Colorize",
                               value = F,
                               status = "primary",
                               shape = "curve",
                               outline = TRUE)
              )
            )
          ),
          selectizeInput(NS(id,"gen_coexp"),
                         "Genes",
                         choices = NULL, 
                         options = list(maxItems = 2,
                                        maxOptions = 20,
                                        placeholder = 'Select two genes for co-expression'),
                         multiple=T)
        ),
        conditionalPanel("typeof output.plot_expression !== 'undefined'", ns = NS(id), 
          box(title = " Co-detection Summary",
              width = NULL,solidHeader = F, collapsible = F, align="center",
            tableOutput(NS(id,"DTCoExp")) %>% withLoader(type='html',loader = 'loader6')
          )
        )
      ),
      column(9,
        box(title = "Scatter Plot",
            width = NULL, solidHeader = T, collapsible = T,
          conditionalPanel("!input.button && typeof output.plot_expression !== 'undefined'", ns = NS(id),
            fluidRow(column=12,align = "left",style='padding-left:12px; padding-right:12px;',
              dropdownButton(
                plotOutput(NS(id,"plot_gradient"),width = "400px",height =  "400px") %>% withSpinner(),
                            circle = TRUE, status = "danger", icon = icon("palette"), width = "300px",
                            tooltip = tooltipOptions(title = "Click to see legend plot",
                                                     placement= "right"),
                            right = F
              )
            )
          ),
          tabsetPanel(id = NS(id,"switcher3"),
                      type = "hidden",
                      selected = "expression_panel",
                      tabPanelBody("cluster_panel",
                                   plotlyOutput(NS(id,"plot_cluster"),height = "100vh") %>% withLoader(type='html',loader = 'dnaspin')
                      ),
                      tabPanelBody("expression_panel",
                                   plotlyOutput(NS(id,"plot_expression"),height = "100vh") %>% withLoader(type='html',loader = 'dnaspin'),
                                   conditionalPanel("typeof output.plot_expression !== 'undefined'", ns = NS(id),
                                   box(title="Co-detection Matrix",width = NULL,
                                       solidHeader = T, collapsible=T,
                                       dropdownButton(
                                         numericInput(NS(id,"pdf_widht_corrplot"),"Widht",value = 7),
                                         numericInput(NS(id,"pdf_heigth_corrplot"),"Heigth",value = 7),
                                         downloadButton(NS(id,'export_corrplot')),
                                         circle = FALSE,
                                         status = "primary",
                                         icon = icon("cog"),
                                         width = "300px",
                                         size= "sm",
                                         up = T,
                                         tooltip = tooltipOptions(title = "Press to Download")
                                       ),
                                       plotOutput(NS(id,"CorrPlot")) %>% withSpinner()
                                   )
                              )
                      )
          )
        ),
      )
    )
  )
}

##### CoExpression Server Module ----
COExpServer <- function(id,sce,point.size=20) {
  moduleServer(id, function(input,output,session) {
    
    ### Observe Events ----
    updatePickerInput(session, 'partitionType', choices = names(colData(sce))[sapply(colData(sce), is.factor)])
    
    updateSelectizeInput(session, 'gen_coexp', choices = rownames(sce), server = TRUE)
    
    dimVector <- reactive({
      sapply(reducedDims(x = sce),FUN = ncol) 
    })
    
    observeEvent(c(dimVector()),{
      updatePickerInput(session,inputId = "DimType", choices = (c("3","2")[c(3,2) %in% dimVector()]))
    })
    
    observeEvent(c(input$DimType,dimVector()), {
      req(!is.null(dimVector()))
      req(input$DimType)
	  updatePickerInput(session,inputId = "plotType", choices = rev(names(which(dimVector() == as.numeric(input$DimType)  | dimVector() > 3))))
    })
    
    observeEvent(ignoreInit = T,input$button,{
      if(input$button) {
        updateTabsetPanel(inputId = "switcher3", selected = "cluster_panel")
      } else{
        updateTabsetPanel(inputId = "switcher3", selected = "expression_panel")
      }
    })
    ### Table and Gen selected ----
    
    CoExpressionL <- reactive({
      req(length(input$gen_coexp) >1) #Two genes min
      genes <- input$gen_coexp
      list(Colors = CoExp_Col(sce,genes = genes),
           Genes = genes)
    })
    
    CoExpressionVtor <- reactive({
      req(length(input$gen_coexp) >1) #Two genes min
      COexp_Vtor(sce = sce,genes = input$gen_coexp)
    })
    
    output$DTCoExp <- renderTable({
      req(!is.null(CoExpressionVtor()))
      data.frame(Genes = CoExpressionVtor()) %>% group_by(Genes) %>% summarise("#" = n(), "%" = (n()/nrow(.))*100)
    })
    
    
    ### Plots ----
    
    OrderPartReact <- eventReactive(input$partitionType,{
      req(input$partitionType)
      Col.and.Order(partition = input$partitionType, sce=sce)
    })
    
        #### Scatter ----
    ClusterPlot <- eventReactive(c(input$plotType,input$partitionType),{
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
                      # customdata= ~colData(sce)[,input$partitionType],
                      size = I(25),span=I(0),text=~colData(sce)[,input$partitionType],hoverinfo='text') %>% 
          toWebGL()
      } else { #2D
        plot_ly(type = "scatter", mode = "markers",source="PlotMix")  %>%
          layout(dragmode = "select",
                 xaxis = list(title = 'Dim1',zeroline=F),
                 yaxis = list(title = 'Dim2',zeroline=F),
                 title = paste(input$partitionType, 'Partition')) %>%
          add_markers(x = ~reducedDim(sce,input$plotType)[,1], y=~reducedDim(sce,input$plotType)[,2], 
                      color = ~I(as.character(OrderPartReact()$colPart[colData(sce)[,input$partitionType]])),
                      name = ~colData(sce)[,input$partitionType],
                      # customdata= ~colData(sce)[,input$partitionType], 
                      size = I(25),span=I(0),text=~colData(sce)[,input$partitionType],hoverinfo='text') %>% 
          toWebGL()
      }
      
    })
    
    
    
    CoExpressionPlot <- eventReactive(c(input$plotType,CoExpressionL(),input$partitionType),{
      req(!is.null(CoExpressionL()))
      req(input$plotType)
      #3D
      if(input$DimType == "3"){
        plot_ly(type = "scatter3d", mode = "markers")  %>%
          layout(dragmode = "select",
                 scene = list(xaxis = list(title = 'Dim1',showgrid=F,visible=F),
                              yaxis = list(title = 'Dim2',showgrid=F,visible=F),
                              zaxis = list(title = 'Dim3',showgrid=F,visible=F)),
                 showlegend = FALSE,
                 title = paste(paste(input$gen_coexp,collapse = " - "), 'Co-expression'),
                 margin = list(l = 0,
                               r = 10,
                               b = 0,
                               t = 40,
                               pad = 0)) %>% 
          add_trace(x = ~reducedDim(sce,input$plotType)[,1], y=~reducedDim(sce,input$plotType)[,2],
                    z=~reducedDim(sce,input$plotType)[,3],
                    size = I(point.size),span=I(0),
                    # name = ~colData(sce)[,"clustersOK"],
                    text= ~colData(sce)[,input$partitionType],
                    hoverinfo = 'text',
                    marker = list(color = ~I(CoExpressionL()$Colors))
                ) %>% 
          toWebGL()
        
      } else { #2D
        plot_ly(type = "scatter", mode = "markers")  %>%
          layout(dragmode = "select",
                 xaxis = list(title = 'Dim1',zeroline=F),
                 yaxis = list(title = 'Dim2',zeroline=F),
                 title = paste(paste(input$gen_coexp,collapse = "-"), 'Co - Localization'),
                 showlegend = FALSE) %>% 
          add_trace(x = ~reducedDim(sce,input$plotType)[,1], y=~reducedDim(sce,input$plotType)[,2],
                    # name = ~colData(sce)[,input$partitionType],
                    text= ~colData(sce)[,input$partitionType],
                    hoverinfo = 'text',
                    size = I(point.size),span=I(0),
                    marker = list(color =I(CoExpressionL()$Colors))) %>% 
          toWebGL()
      }
      
    })
    
    output$plot_cluster <- renderPlotly({
      req(!is.null(ClusterPlot()))
        ClusterPlot()
    })
    
    output$plot_expression <- renderPlotly({
      req(!is.null(CoExpressionPlot()))
      CoExpressionPlot()
    })
    
     output$plot_gradient <- renderPlot({
       req(length(input$gen_coexp) >1)
       # req(input$button == F)
       plot2dgradient(gen1 = input$gen_coexp[1],gen2 = input$gen_coexp[2])
     })
    
      
        #### CorrPlot ----
     
     
     CorrPlot <- reactive({
       req(!is.null(CoExpressionVtor()))
       req(input$partitionType)
       apply(table(CoExpressionVtor(),colData(sce)[,input$partitionType]),2,function(x){(x/sum(x))}) %>% 
         .[nrow(.):1,] %>% #Reverse the matrix to make it similar to Table 
         t %>% 
         ggcorrplot::ggcorrplot(lab=T) +
         scale_fill_gradient2(midpoint = 0.5,low="white",high ="#E46726",mid="#ffa18c") + labs(fill="Proportion")
     })
     
     output$CorrPlot <- renderPlot({
       req(!is.null(CorrPlot()))
       CorrPlot()
     })
     
     output$export_corrplot = downloadHandler(
       filename = function() {"CoDetection_Matrix.pdf"},
       content = function(file) {
         pdf(file,
             width = input$pdf_widht_corrplot,
             height = input$pdf_heigth_corrplot
         )
         CorrPlot() %>% plot()
         dev.off()
       }
     )
     
  }) 
}
