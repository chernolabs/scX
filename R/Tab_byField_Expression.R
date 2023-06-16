#####                              ######## 
#####         Expression Tab       ######## 
#####                              ######## 

##### Expression UI Module ----
Numeric_ExpressionUI <- function(id) {
  tagList(
    fluidRow(
      column(3,
        box(title = htmltools::span(icon("gears"), " Settings"),
            width = NULL, status = "primary",solidHeader = T,collapsible = F,
            fluidRow(
              column(12,style='padding-left:12px; padding-right:12px;',
                    align="center",
                     pickerInput(inputId = NS(id,"numericType"), 
                                 label = "Numeric",
                                 choices = NULL)
              )
          ),
          conditionalPanel("input.scatter_heatmap == 'scatter'",ns=NS(id),
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
          conditionalPanel("(typeof output.plot_expression !== 'undefined' || input.scatter_heatmap == 'heatmap') && input.scatter_heatmap !== 'MultiLines'", ns = NS(id),
            fluidRow(
              column(8,style='padding-left:12px; padding-right:3px;',
                pickerInput(inputId = NS(id,"partitionType"), 
                            label = "Partition",
                            choices = NULL)
              ),
              conditionalPanel("input.scatter_heatmap == 'scatter'",ns=NS(id),
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
            )
          ),
          conditionalPanel("input.scatter_heatmap == 'heatmap'",ns=NS(id),
            prettySwitch(NS(id,"cluster_row"),
                         "Cluster Row",
                         value = F,
                         status = "primary",
                         fill = TRUE),
			      prettySwitch(NS(id,"norm_heat"),
                         "Norm per gene",
                         value = F,
                         status = "primary",
                         fill = TRUE),
			      conditionalPanel("input.partitionType != 'None'",ns=NS(id),
			         prettySwitch(NS(id,"split_column"),
			                      "Split by Column",
			                      value = F,
			                      status = "primary",
			                      fill = TRUE)
			      )
          ),
          # conditionalPanel("input.scatter_heatmap == 'MultiLines'",ns=NS(id),
          #   prettySwitch(NS(id,"ord_MultiLines"),
          #                "Cluster Row",
          #                value = F,
          #                status = "primary",
          #                fill = TRUE)
          # ),
          fluidRow(
            column(12,style='padding-left:12px; padding-right:12px;',align="center",
              switchInput(NS(id,"GL_T"), 
                          label = "Upload GeneList",
                          size = "small",
                          width = NULL,
                          labelWidth = "100px"),
            )
          ),
          tabsetPanel(id = NS(id,"switcher"),
                      type = "hidden",
                      selected = "panel1",
            tabPanelBody("panel1",
              fluidRow(
                column(11,style='padding-left:0px; padding-right:2px;',
                  selectizeInput(NS(id,"gen_exp"),
                                 label=NULL,
                                 choices = NULL, 
                                 options = list(maxItems = 10,
                                                maxOptions = 20,
                                                placeholder = 'Please select genes to plot'),
                                 width = NULL,
                                 multiple=T)
                ),
                column(1,style='padding-left:2px; padding-right:2px; padding-top:4px',
                  actionBttn((NS(id,"action")),
                             label = NULL,
                             style = "unite",
                             color = "primary",
                             size = "xs",
                             icon = icon("play"))
                )
              )
            ),
            tabPanelBody("panel2",
              fluidRow(
                column(12,style='padding-left:0px; padding-right:0px;',
                  fileInput(NS(id,'listGenes'),
                            label = NULL,
                            multiple = T, 
                            accept = c("txt/csv", "text/comma-separated-values,
                                        text/plain", ".csv", ".xlsx"),
                            buttonLabel = "Search",
                            placeholder = "Select Gene list"),
                  htmlOutput(NS(id,"missingGenes"))
                )
              )
            )
          )
        )
      ),
      column(9,
        tabBox(id = NS(id,"scatter_heatmap"),
               selected = "scatter",
               width = NULL,
          tabPanel("Scatter",value = "scatter",
            box(title = "Scatter Plot",
                width = NULL, solidHeader = T, collapsible = T,
                footer = tagList(shiny::icon("cat"), "Nya"),
                tabsetPanel(id = NS(id,"switcher3"),
                            type = "hidden",
                            selected = "expression_panel",
                            tabPanelBody("cluster_panel",
                                         plotlyOutput(NS(id,"plot_cluster"),height = "100vh") %>% withLoader(type='html',loader = 'dnaspin')
                            ),
                            tabPanelBody("expression_panel",
                                         plotlyOutput(NS(id,"plot_expression"),height = "100vh") %>% withLoader(type='html',loader = 'dnaspin'),
                                         uiOutput(NS(id,"Lines.Bar_Input")),
                                         conditionalPanel("typeof output.plot_expression !== 'undefined'", ns = NS(id),
                                         tabsetPanel(id = NS(id,"switcher2"),
                                                     type = "hidden",
                                                     selected = "Lines_panel",
                                                     tabPanelBody("Lines_panel",
                                                                  dropdownButton(
																	fluidRow(
																		column(4, style='padding-left:12px; padding-right:3px;', numericInput(NS(id,"pdf_width_linesplot"),"Width",value = 7)),
																		column(4, style='padding-left:3px; padding-right:3px;', numericInput(NS(id,"pdf_height_linesplot"),"Height",value = 7)),
																		column(4, style='padding-left:3px; padding-right:12px; padding:16px', downloadButton(NS(id,'export_linesplot')))
																	),
                                                                    circle = FALSE,
                                                                    status = "primary",
                                                                    icon = icon("fa-thin fa-download"),
                                                                    width = "300px",
                                                                    size= "sm",
                                                                    up = T,
                                                                    tooltip = tooltipOptions(title = "Press to Download")
                                                                  ),
                                                                  plotlyOutput(NS(id,"plot_Lines")) %>% withSpinner(),
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
                                                                    icon = icon("fa-thin fa-download"),
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
          ),
          tabPanel("Heatmap", value= "heatmap",
            box(width = NULL, solidHeader = T, collapsible = F,
                footer = tagList(shiny::icon("cat"), "Nya"),
                dropdownButton(
					fluidRow(
						column(4, style='padding-left:12px; padding-right:3px;', numericInput(NS(id,"pdf_width_heatmap"),"Width",value = 7)),
						column(4, style='padding-left:3px; padding-right:3px;', numericInput(NS(id,"pdf_height_heatmap"),"Height",value = 7)),
						column(4, style='padding-left:3px; padding-right:12px; padding:16px', downloadButton(NS(id,'export_heatmap')))
					),
                  circle = FALSE,
                  status = "primary",
                  icon = icon("fa-thin fa-download"),
                  width = "300px",
                  size= "sm",
                  up = F,
                  tooltip = tooltipOptions(title = "Press to Download")
                ),
              plotOutput(NS(id,"plot_heatmap"),height = "100vh") %>% withLoader(type='html',loader = 'dnaspin')
            )
          ),
          tabPanel("MultiLines", value= "MultiLines",
            box(width = NULL,solidHeader = T,collapsible = F,
                footer = tagList(shiny::icon("cat"), "Nya"),
                dropdownButton(
					fluidRow(
						column(4, style='padding-left:12px; padding-right:3px;', numericInput(NS(id,"pdf_width_multilines"),"Width",value = 7)),
						column(4, style='padding-left:3px; padding-right:3px;', numericInput(NS(id,"pdf_height_multilines"),"Height",value = 7)),
						column(4, style='padding-left:3px; padding-right:12px; padding:16px', downloadButton(NS(id,'export_multilines')))
					),
                  circle = FALSE,
                  status = "primary",
                  icon = icon("fa-thin fa-download"),
                  width = "300px",
                  size= "sm",
                  up = F,
                  tooltip = tooltipOptions(title = "Press to Download")
                ),
              plotlyOutput(NS(id,"plot_MultiLines"),height = "100vh") %>% withLoader(type='html',loader = 'dnaspin')
            )
          )
        )
      )
    )
  )
}


##### Expression Server Module ----
Numeric_ExpressionServer <- function(id,sce,point.size=20) {
  moduleServer(id, function(input,output,session) {
    ### Observe Events ----
    observeEvent(ignoreInit = T,input$GL_T,{
      if(input$GL_T) {
        updateTabsetPanel(inputId = "switcher", selected = "panel2")
      } else{
        updateTabsetPanel(inputId = "switcher", selected = "panel1")
      }
    })
    
    observeEvent(ignoreInit = T,input$button,{
      if(input$button) {
        updateTabsetPanel(inputId = "switcher3", selected = "cluster_panel")
      } else{
        updateTabsetPanel(inputId = "switcher3", selected = "expression_panel")
      }
    })
    
    updatePickerInput(session, 'partitionType', 
                         choices = c('None',names(colData(sce))[sapply(colData(sce), is.factor)]))
    
    updatePickerInput(session, 'numericType', 
                      choices = names(colData(sce))[sapply(colData(sce), is.numeric)])
    
    updateSelectizeInput(session, 'gen_exp', choices = rownames(sce), server = TRUE)
    
    dimVector <- reactive({
        sapply(reducedDims(x = sce),FUN = ncol) 
    })
    
    observeEvent(c(dimVector()),{
      #req(input$scatter_heatmap == "scatter")
      updatePickerInput(session,inputId = "DimType", choices = (c("3","2")[c(3,2) %in% dimVector()]))
    })
    
    observeEvent(c(input$DimType,dimVector()), {
      #req(input$scatter_heatmap == "scatter")
      req(!is.null(dimVector()))
      req(input$DimType)
	  updatePickerInput(session,inputId = "plotType", choices = rev(names(which(dimVector() == as.numeric(input$DimType)  | dimVector() > 3))))
    })
    
    observeEvent(input$Cell_Exp,{
      req(input$Cell_Exp)
      switch(input$Cell_Exp,
             'Lines'    = updateTabsetPanel(inputId = "switcher2",  selected = "Lines_panel"),
             'SpikePlot' = updateTabsetPanel(inputId = "switcher2", selected = "SpikePlot_panel")
      )
    })
    ### Gen selected & data preparation ----
        #### Selected Gene -----
    genes.L <- eventReactive(input$action,{
      req(input$GL_T == F)
      req(input$gen_exp)
      input$gen_exp
    })
    
    ExpressionL <- eventReactive(genes.L(),{
      #req(input$scatter_heatmap == "scatter")
      req(input$GL_T == F)
      req(!is.null(genes.L()))
      if(length(genes.L()) > 1) {
        exp_vtor <- apply(logcounts(sce)[genes.L(),],2,mean)
      } else {
        exp_vtor <- logcounts(sce)[genes.L(),]
      }
      list(Exp = exp_vtor,Genes=genes.L()) #these because i dont want to change the gene.L object.
    })
    
    HeatmapL <- eventReactive(c(genes.L(),input$norm_heat),{
      #req(input$scatter_heatmap == "heatmap")
      req(input$GL_T == F)
      req(!is.null(genes.L()))
      value <- ifelse(input$norm_heat,yes = "logcounts.norm",no = "logcounts") #To swtich between norm ot not normalize expression gene.
      exp_mtx <-as.matrix(assay(sce,value)[genes.L(),])
      if(length(genes.L()) == 1) {
        exp_mtx <- t(exp_mtx)
      }
      rownames(exp_mtx) <- genes.L()
      
      exp_mtx
    })
    
    #### Gene List Uploaded  -----
    genes.GL <- eventReactive(input$listGenes,{
      req(input$listGenes)
      GL <- genesList(dataPath = input$listGenes)
      no.genes <- GL[!(GL %in% rownames(sce))] #Keep the missing values
      genes <- GL[(GL %in% rownames(sce))]
      
      list(genes = genes,miss= no.genes)
    })
    
    output$missingGenes <- renderText({
      req(length(genes.GL()$miss)>0)
      paste('<b style="color:red;">The following genes were not found in the dataset:<b><br>', 
            paste(genes.GL()$miss,
                  collapse = ", ")
      )
    })
    
    ExpressionGL <- eventReactive(genes.GL(),{
      #req(input$scatter_heatmap == "scatter")
      req(input$GL_T == T)
      req(length(genes.GL()$genes) > 0)
      
      if(length(genes.GL()$genes) > 1) {
        exp_vtor <- apply(logcounts(sce)[genes.GL()$genes,],2,mean)
      } else {
        exp_vtor <- logcounts(sce)[genes.GL()$genes,]
      }
      
      list(Exp = exp_vtor,Genes=genes.GL()$genes)
    })
    
    HeatmapGL <- eventReactive(c(genes.GL(),input$norm_heat),{
      #req(input$scatter_heatmap == "heatmap")
      req(input$GL_T == T)
      req(length(genes.GL()$genes) > 0)
      value <- ifelse(input$norm_heat,yes = "logcounts.norm",no = "logcounts") #To swtich between norm ot not normalize expression gene.
      exp_mtx <-as.matrix(assay(sce,value)[genes.GL()$genes,])
      if(length(genes.GL()$genes) == 1) {
        exp_mtx <- t(exp_mtx)
      }
      rownames(exp_mtx) <- genes.GL()$genes
      
      exp_mtx
    })
    
    #### Selecting ----
    #I do that because if it change, it haven't to calculated everything again, just change the dataframe
    HeatmapF <- reactive({
      #req(input$scatter_heatmap == "heatmap")
      if(input$GL_T){
        HeatmapGL()    
      } else { 
        HeatmapL()
      }
    })
    
    ExpressionF <- reactive({
      if(input$GL_T){
        ExpressionGL()    
      } else { 
        ExpressionL()
      }
    })
    
    
    #### Plots ----
    
    OrderPartReact <- eventReactive(input$partitionType,{
      req(input$partitionType != 'None')
      Col.and.Order(partition = input$partitionType, sce=sce)
    })
    
    ### Scatter ----
    ClusterPlot <- eventReactive(c(input$plotType,input$partitionType,input$numericType),{
      #req(input$scatter_heatmap == "scatter")
      #3D
      if(input$DimType == "3"){
        plot_ly(type = "scatter3d", mode = "markers",source = "PlotMix",colors = 'YlOrRd')  %>%
          layout(dragmode = "select",
                 scene = list(xaxis = list(title = 'Dim1',showgrid=F,visible=F),
                              yaxis = list(title = 'Dim2',showgrid=F,visible=F),
                              zaxis = list(title = 'Dim3',showgrid=F,visible=F)),
                 title = input$numeircType,
                 legend= list(x=1,y=1),
                 showlegend = T, #(input$partitionType != 'None'),
                 margin = list(l = 0,
                               r = 10,
                               b = 0,
                               t = 40,
                               pad = 0)) %>% 
          add_markers(x=~reducedDim(sce,input$plotType)[,1],
                      y=~reducedDim(sce,input$plotType)[,2],
                      z=~reducedDim(sce,input$plotType)[,3],
                      text= if(input$partitionType == 'None') {NULL} else {~colData(sce)[,input$partitionType]},
                      hoverinfo = 'text',
                      color = ~colData(sce)[,input$numericType],
                      name = if(input$partitionType == 'None') {NULL} else {~colData(sce)[,input$partitionType]},
                      size = I(point.size),
                      span=I(0)) %>% 
          colorbar(title =input$numericType,x=0,y=1) %>% 
          toWebGL()
      } else { #2D
        plot_ly(type = "scatter", mode = "markers",source="PlotMix",colors = 'YlOrRd')  %>%
          layout(dragmode = "select",
                 xaxis = list(title = 'Dim1',zeroline=F),
                 yaxis = list(title = 'Dim2',zeroline=F),
                 showlegend = F,
                 title = input$numericType) %>%
          add_markers(x = ~reducedDim(sce,input$plotType)[,1],
                      y=~reducedDim(sce,input$plotType)[,2], 
                      color = ~colData(sce)[,input$numericType],
                      text=  if(input$partitionType == 'None') {NULL} else {~colData(sce)[,input$partitionType]},
                      hoverinfo = 'text',
                      size = I(point.size),
                      span=I(0)) %>% 
          colorbar(title = input$numericType) %>% 
          toWebGL()
      }
      
    })
    
    ExpressionPlot <- eventReactive(c(input$plotType,ExpressionF(),input$partitionType),{
      req(!is.null(ExpressionF()))
      req(input$plotType)
      #3D
      if(input$DimType == "3"){
        plot_ly(type = "scatter3d", mode = "markers")  %>%
          layout(dragmode = "select",
                 scene = list(xaxis = list(title = 'Dim1',showgrid=F,visible=F),
                              yaxis = list(title = 'Dim2',showgrid=F,visible=F),
                              zaxis = list(title = 'Dim3',showgrid=F,visible=F)),
                 legend= list(x=1,y=1),
                 showlegend = T, #(input$partitionType != 'None'),
                 title = ifelse(length(ExpressionF()$Genes)>1,'Mean Expression', 'Expression'),
                 margin = list(l = 0,
                               r = 10,
                               b = 0,
                               t = 40,
                               pad = 0)) %>% 
          add_markers(x = ~reducedDim(sce,input$plotType)[,1], y=~reducedDim(sce,input$plotType)[,2],
                      z=~reducedDim(sce,input$plotType)[,3],
                      text= if(input$partitionType == 'None') {NULL} else {~colData(sce)[,input$partitionType]},
                      hoverinfo = 'text',
                      color = ~ExpressionF()$Exp,
                      name = if(input$partitionType == 'None') {NULL} else {~colData(sce)[,input$partitionType]},
                      size = I(point.size),span=I(0)) %>% 
          colorbar(title = "log(counts)",x=0,y=1) %>% 
          toWebGL()
        
      } else { #2D
        plot_ly(type = "scatter", mode = "markers")  %>%
          layout(dragmode = "select",
                 xaxis = list(title = 'Dim1',zeroline=F),
                 yaxis = list(title = 'Dim2',zeroline=F),
                 title = ifelse(length(ExpressionF()$Genes)>1,'Mean Expression', 'Expression'),
                 showlegend = FALSE) %>% 
          add_markers(x = ~reducedDim(sce,input$plotType)[,1], y=~reducedDim(sce,input$plotType)[,2],
                      color = ~ExpressionF()$Exp, 
                      size = I(point.size),span=I(0),
                      text=  if(input$partitionType == 'None') {NULL} else {~colData(sce)[,input$partitionType]},
                      #name = ~colData(sce)[,input$partitionType],
                      hoverinfo = 'text') %>% 
          colorbar(title = "log(counts)") %>% 
          toWebGL()
      }
      
    })
    
    output$plot_cluster <- renderPlotly({
      req(!is.null(ClusterPlot()))
      ClusterPlot()
    })
    
    output$plot_expression <- renderPlotly({
      req(!is.null(ExpressionPlot()))
      ExpressionPlot()
    })
    
    
    ### Heatmaps ----
    #Here I made everything in the same reactive object to manipulated all the reactive order at the same time
    Heatmap_Plot <- eventReactive(c(HeatmapF(),input$partitionType,input$numericType,input$cluster_row,input$split_column),{
      #req(input$scatter_heatmap == "heatmap")
      req(input$partitionType)
      req(input$numericType)
      req(!is.null(HeatmapF()))
      # req(!is.null(OrderPartReact()))
      
      col_fun = circlize::colorRamp2(c(0, max(colData(sce)[,input$numericType,drop=T], na.rm = T)), hcl_palette ='YlOrRd',reverse = T)
      ht <-   HeatmapAnnotation(numeric = colData(sce)[,input$numericType],
                                Type = if(input$partitionType == 'None') {NULL} else {colData(sce)[,input$partitionType]},
                                col= if(input$partitionType == 'None') {list(numeric = col_fun)} else {list(numeric = col_fun,Type=OrderPartReact()$colPart)} ,
                                show_legend = c(Type = FALSE),
                                annotation_label = if(input$partitionType == 'None') {c(input$numericType)} else {c(input$numericType,input$partitionType)},
                                show_annotation_name = T)
      
      split_col <- NULL
      if(input$partitionType != 'None' & input$split_column){
        split_col <- colData(sce)[,input$partitionType]
        }
        h1 <- Heatmap(as.matrix(HeatmapF()),
                      col = if(max(HeatmapF())==min(HeatmapF())) {viridis::viridis(1)} else {viridis::viridis(100)},
                      border =F,
                      name = "Gene expression",
                      cluster_rows = input$cluster_row,
                      cluster_columns = F,
                      row_names_side = "left",
                      column_title_rot = 45,
                      column_title_side = "bottom",
                      row_title = "Genes",
                      row_gap = unit(1, "mm"),
                      column_gap = unit(1, "mm"),
                      show_row_names = T,
                      show_column_names = F,
                      top_annotation = ht,
                      column_order = order(colData(sce)[,input$numericType]),
                      column_split = split_col,
                      cluster_column_slices = F,
                      use_raster = TRUE,
                      raster_by_magick = TRUE
        )
      h1
    })

    output$plot_heatmap <- renderPlot({
      #req(input$scatter_heatmap == "heatmap")
      req(!is.null(Heatmap_Plot()))
      Heatmap_Plot()
    })
    
    ####  Lines&SpikePlots ----
     output$Lines.Bar_Input <- renderUI({
       req(!is.null(ExpressionF()))
       radioGroupButtons(inputId = NS(id,"Cell_Exp"), label=NULL,choices = c("Lines", "SpikePlot"),
                         direction = "horizontal",justified = T,individual=T)
     })
     
     LinesReact <- reactive({
       req(!is.null(ExpressionF()))
       
       p <- data.frame(Exp=ExpressionF()$Exp,Numeric = colData(sce)[,input$numericType])
       if(input$partitionType != 'None') {
         p$partitionType <- colData(sce)[,input$partitionType]
       }
       p
     })
     
     LinesPlot  <- reactive({
       #req(input$scatter_heatmap == "scatter")
       req(!is.null(LinesReact()))
        p <- LinesReact() %>% ggplot(aes(y=Exp,x=Numeric)) + geom_point(aes(col = if(input$partitionType != 'None'){partitionType} else{NULL}),alpha=0.5) + geom_smooth(se = F) +
         labs(color = input$partitionType) + 
          xlab(input$numericType) + 
          ylab("log(counts)") +
          ggtitle(ifelse(length(ExpressionF()$Genes)>1,"mean expression","expression"))  + 
          scale_color_manual(values=if(input$partitionType != 'None'){ OrderPartReact()$colPart} else{'grey'})
        p 
     })
     
     output$plot_Lines  <- renderPlotly({
       #req(input$scatter_heatmap == "scatter")
       req(!is.null(LinesPlot()))
       LinesPlot() %>% ggplotly() %>% toWebGL()
     })
     
     SpikePlot <-reactive({
       #req(input$scatter_heatmap == "scatter")
       req(!is.null(ExpressionF()))
       if(input$partitionType == 'None'){
         ord <- order(colData(sce)[,input$numericType])
         barplot(ExpressionF()$Exp[ord],
                      col = "black",
                      border = "black",
                      ylab = "log(counts)",
                 main = ifelse(length(ExpressionF()$Genes)>1,"mean expression","expression"),
                 names.arg = F) 
         # legend("bottom", legend = names(OrderPartReact()$colPart), col = OrderPartReact()$colPart,
                # pch=19, ncol=6, xpd=T, inset=c(0,-0.25))
         abline(h=mean(ExpressionF()$Exp[ExpressionF()$Exp>0]),lty=2,col="grey")
       } else{
         ord <- order(ordered(colData(sce)[,input$partitionType]), colData(sce)[,input$numericType])
         m <- barplot(ExpressionF()$Exp[ord],
                      col = OrderPartReact()$colPart[colData(sce)[,input$partitionType]][ord],
                      border = OrderPartReact()$colPart[colData(sce)[,input$partitionType]][ord],
                      ylab = "log(counts)", main = ifelse(length(ExpressionF()$Genes)>1,"mean expression","expression"), names.arg = F) 
         legend("bottom", legend = names(OrderPartReact()$colPart), col = OrderPartReact()$colPart,
                pch=19, ncol=6, xpd=T, inset=c(0,-0.25))
         lines(x = m,
               tapply(ExpressionF()$Exp,
                      INDEX = colData(sce)[,input$partitionType],
                      FUN = mean)[colData(sce)[,input$partitionType]][ord],
               lty=2,col="black")
       }
       graph <- recordPlot()
       graph
     })
     
     output$plot_SpikePlot <- renderPlot({
       req(!is.null(SpikePlot()))
       SpikePlot() %>% print()
     })
     
     #### Line plot MultiGenes ----
     MultiLinesPlot <- reactive({
       #req(input$scatter_heatmap == "MultiLines")
       req(input$numericType)
       if(input$GL_T){
         req(length(genes.GL()$genes) > 0)
         feature <- genes.GL()$genes
       } else { 
         req(!is.null(genes.L()))
         feature <- genes.L() 
       }
       
       df_plot <- assay(sce,"logcounts")[feature,,drop=F] %>% as.matrix %>% t %>% as.data.frame
       # Add cell ID and identity classes
       df_plot$Cell <- rownames(df_plot)
       df_plot$Numeric <- colData(sce)[,input$numericType]
       df_plot <- reshape2::melt(df_plot, id.vars = c("Cell","Numeric"), measure.vars = feature,
                                variable.name = "Feat", value.name = "Expr")
       
       p <- df_plot %>% ggplot(aes(y=Expr,x=Numeric,col = Feat)) + geom_smooth(se = F) +
         xlab(input$numericType) + 
         ylab("log(counts)") +
         ggtitle("Expression")   
       p 
     })
     
     output$plot_MultiLines <- renderPlotly({
       req(!is.null(MultiLinesPlot()))
       MultiLinesPlot() %>% ggplotly() %>% toWebGL()
     })
     
     ### Downloads -----
     
     output$export_SpikePlot = downloadHandler(
       filename = function() {"SpikePlot_Fields.pdf"},
       content = function(file) {
         pdf(file,
             width = input$pdf_width_SpikePlot,
             height = input$pdf_height_SpikePlot
         )
         SpikePlot() %>% print()
         dev.off()
       })
     
     output$export_heatmap = downloadHandler(
       filename = function() {"Heatmap_Fields.pdf"},
       content = function(file) {
         pdf(file,
             width = input$pdf_width_heatmap,
             height = input$pdf_height_heatmap
         )
         Heatmap_Plot() %>% plot()
         dev.off()
       }
     )
     
     output$export_linesplot = downloadHandler(
       filename = function() {"LinesPlot_Fields.pdf"},
       content = function(file) {
         pdf(file,
             width = input$pdf_width_linesplot,
             height = input$pdf_height_linesplot
         )
         LinesPlot() %>% plot()
         dev.off()
       }
     )
     
     output$export_multilines = downloadHandler(
       filename = function() {"MultiLines_Fields.pdf"},
       content = function(file) {
         pdf(file,
             width = input$pdf_width_multilines,
             height = input$pdf_height_multilines
         )
         MultiLinesPlot() %>% plot()
         dev.off()
       }
     )
  }) 
}

