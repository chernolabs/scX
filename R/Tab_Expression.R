#####                              ######## 
#####         Expression Tab       ######## 
#####                              ######## 

##### Expression UI Module ----
ExpressionUI <- function(id) {
  tagList(
    fluidRow(
      column(3,
        box(title = htmltools::span(icon("gears"), " Settings"),
			width = NULL, status = "primary",solidHeader = T,collapsible = F,
			conditionalPanel("!input.button || input.scatter_heatmap != 'scatter'",ns=NS(id),
			fluidRow( # gene list upload off/on
				column(12,style='padding-left:12px; padding-right:12px;',align="center",
					switchInput(NS(id,"GL_T"), 
						label = "Upload Gene List",
						size = "small",
						width = NULL,
						labelWidth = "100px"),
				)
			),
			# gene selection or list
			tabsetPanel(id = NS(id,"switcher"),
				type = "hidden",
				selected = "panel1",
				tabPanelBody("panel1",
					fluidRow(
						column(11,style='padding-left:0px; padding-right:2px;',
							selectizeInput(NS(id,"gen_exp"),
							label=NULL,
							choices = NULL, 
							options = list(maxOptions = 20,
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
								text/plain", ".csv", ".xlsx", ".xls"),
							buttonLabel = "Search",
							placeholder = "Select Gene list"),
							htmlOutput(NS(id,"missingGenes"))
						)
					)
				)
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
          #conditionalPanel("typeof output.plot_expression !== 'undefined' || input.scatter_heatmap == 'heatmap' || input.scatter_heatmap == 'dotplot' || input.scatter_heatmap == 'stackVln'", ns = NS(id),
            fluidRow(
              column(8,style='padding-left:12px; padding-right:3px;',
                pickerInput(inputId = NS(id,"partitionType"), 
                            label = "Category",
                            choices = NULL)
              ),
              conditionalPanel("input.scatter_heatmap == 'scatter'",ns=NS(id),
                column(4,style='padding-left:3px; padding-right:1px;padding-top:12px',
                  br(),
                  prettyCheckbox(NS(id,"button"),
                                 label="Show",
                                 value = F,
                                 status = "primary",
                                 shape = "curve",
                                 outline = TRUE)
                )
              )
            ),
          #),
          conditionalPanel("input.scatter_heatmap == 'heatmap'",ns=NS(id),
            prettySwitch(NS(id,"cluster_row"),
                         "Cluster Row",
                         value = F,
                         status = "primary",
                         fill = TRUE),
            prettySwitch(NS(id,"cluster_column"),
                         "Cluster Column",
                         value = F,
                         status = "primary",
                         fill = TRUE),
			prettySwitch(NS(id,"norm_heat"),
                         "Norm per gene",
                         value = F,
                         status = "primary",
                         fill = TRUE),
            prettySwitch(NS(id,"mean_heat"),
                         "Mean per group",
                         value = T,
                         status = "primary",
                         fill = TRUE)
            
          ),
          conditionalPanel("input.scatter_heatmap == 'dotplot'",ns=NS(id),
            prettySwitch(NS(id,"ord_dotplot"),
                         "Cluster Row",
                         value = F,
                         status = "primary",
                         fill = TRUE),
            prettySwitch(NS(id,"scale_dotplot"),
                         "Scale",
                         value = F,
                         status = "primary",
                         fill = TRUE),
            prettySwitch(NS(id,"center_dotplot"),
                         "Center",
                         value = F,
                         status = "primary",
                         fill = TRUE)
          ),
          conditionalPanel("input.scatter_heatmap == 'stackVln'",ns=NS(id),
            prettySwitch(NS(id,"ord_stackVln"),
                         "Cluster Row",
                         value = F,
                         status = "primary",
                         fill = TRUE)
          )
		)
      ),
      column(9,
        tabBox(id = NS(id,"scatter_heatmap"),
               selected = "scatter",
               width = NULL,
          tabPanel("Scatter",value = "scatter",
            box(#title = "Scatter Plot",
                width = NULL, solidHeader = T, collapsible = F,
                footer = tagList(shiny::icon("cat"), "Nya"),
                tabsetPanel(id = NS(id,"switcher3"),
                            type = "hidden",
                            selected = "expression_panel",
                tabPanelBody("cluster_panel",
                  plotlyOutput(NS(id,"plot_cluster"),height = "80vh") %>% withLoader(type='html',loader = 'dnaspin')
                ),
                tabPanelBody("expression_panel",
                             plotlyOutput(NS(id,"plot_expression"),height = "80vh") %>% withLoader(type='html',loader = 'dnaspin'),
                             uiOutput(NS(id,"Violin.Bar_Input")),
                             conditionalPanel("typeof output.plot_expression !== 'undefined'", ns = NS(id),
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
                  icon = icon("download"),
                  width = "300px",
                  size= "sm",
                  up = F,
                  tooltip = tooltipOptions(title = "Download")
                ),
              plotOutput(NS(id,"plot_heatmap"),height = "80vh") %>% withLoader(type='html',loader = 'dnaspin')
            )
          ),
          tabPanel("DotPlot", value= "dotplot",
            box(width = NULL,solidHeader = T,collapsible = F,
                footer = tagList(shiny::icon("cat"), "Nya"),
                dropdownButton(
                  fluidRow(
						column(4, style='padding-left:12px; padding-right:3px;', numericInput(NS(id,"pdf_width_dotplot"),"Width",value = 7)),
						column(4, style='padding-left:3px; padding-right:3px;', numericInput(NS(id,"pdf_height_dotplot"),"Height",value = 7)),
						column(4, style='padding-left:3px; padding-right:12px; padding:16px', downloadButton(NS(id,'export_dotplot')))
					),
                  circle = FALSE,
                  status = "primary",
                  icon = icon("download"),
                  width = "300px",
                  size= "sm",
                  up = F,
                  tooltip = tooltipOptions(title = "Download")
                ),
              plotlyOutput(NS(id,"plot_DotPlot"),height = "80vh") %>% withLoader(type='html',loader = 'dnaspin')
            )
          ),
          tabPanel("StackedViolin", value= "stackVln",
            box(width = NULL,solidHeader = T,collapsible = F,
                footer = tagList(shiny::icon("cat"), "Nya"),
                dropdownButton(
					fluidRow(
						column(4, style='padding-left:12px; padding-right:3px;', numericInput(NS(id,"pdf_width_StackedViolin"),"Width",value = 7)),
						column(4, style='padding-left:3px; padding-right:3px;', numericInput(NS(id,"pdf_height_StackedViolin"),"Height",value = 7)),
						column(4, style='padding-left:3px; padding-right:12px; padding:16px', downloadButton(NS(id,'export_StackedViolin')))
					),
                  circle = FALSE,
                  status = "primary",
                  icon = icon("download"),
                  width = "300px",
                  size= "sm",
                  up = F,
                  tooltip = tooltipOptions(title = "Download")
                ),
              plotOutput(NS(id,"plot_stackVln"),height = "80vh") %>% withLoader(type='html',loader = 'dnaspin')
            )
          )
        )
      )
    )
  )
}


##### Expression Server Module ----
ExpressionServer <- function(id,sce,point.size=20) {
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
                         choices = names(colData(sce))[sapply(colData(sce), is.factor)])
    
    updateSelectizeInput(session, 'gen_exp', choices = rownames(sce), server = TRUE)
    
    dimVector <- reactive({
        sapply(reducedDims(x = sce),FUN = ncol) 
    })
    
    observeEvent(dimVector(),{
      if(any(dimVector() > 3)){
        opt <- c('3','2')
      } else{
        opt <- (c("3","2")[c(3,2) %in% dimVector()])
      }
      updatePickerInput(session,inputId = "DimType", choices =opt)
    })
    
    observeEvent(input$DimType, {
      req(input$DimType)
	  updatePickerInput(session,inputId = "plotType", choices = rev(names(which(dimVector() == as.numeric(input$DimType)  | dimVector() > 3))))
    })
    
    observeEvent(input$Cell_Exp,{
      switch(input$Cell_Exp,
             'Violin'    = updateTabsetPanel(inputId = "switcher2",  selected = "Violin_panel"),
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
      req(input$GL_T == F)
      req(!is.null(genes.L()))
      if(length(genes.L()) > 1) {
        exp_vtor <- colSums(logcounts(sce)[genes.L(),])/length(genes.L())
      } else {
        exp_vtor <- logcounts(sce)[genes.L(),]
      }
      list(Exp = exp_vtor,Genes=genes.L()) #these because i dont want to change the gene.L object.
    })
    
    HeatmapL <- eventReactive(c(genes.L(),input$norm_heat),{
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
      req(input$GL_T == T)
      req(length(genes.GL()$genes) > 0)
      
      if(length(genes.GL()$genes) > 1) {
        exp_vtor <- colSums(logcounts(sce)[genes.GL()$genes,])/length(genes.GL()$genes)
      } else {
        exp_vtor <- logcounts(sce)[genes.GL()$genes,]
      }
      
      list(Exp = exp_vtor,Genes=genes.GL()$genes)
    })
    
    HeatmapGL <- eventReactive(c(genes.GL(),input$norm_heat),{
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
      req(input$partitionType)
      Col.and.Order(partition = input$partitionType, sce=sce)
    })
    
    ### Scatter ----
    output$plot_cluster <- renderPlotly({
      #3D
      if(isolate(input$DimType) == "3"){
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
                      size = I(point.size),span=I(0),text=~colData(sce)[,input$partitionType],hoverinfo='text') %>% 
          suppressWarnings(toWebGL())
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
                      size = I(point.size),span=I(0),text=~colData(sce)[,input$partitionType],hoverinfo='text') %>% 
          toWebGL()
      }
      
    })
    
    output$plot_expression <- renderPlotly({
      req(!is.null(ExpressionF()))
      req(input$plotType)
      #3D
      if(isolate(input$DimType) == "3"){
        plot_ly(type = "scatter3d", mode = "markers")  %>%
          layout(dragmode = "select",
                 scene = list(xaxis = list(title = 'Dim1',showgrid=F,visible=F),
                              yaxis = list(title = 'Dim2',showgrid=F,visible=F),
                              zaxis = list(title = 'Dim3',showgrid=F,visible=F)),
                 legend= list(x=1,y=1),
                 showlegend = TRUE,
                 title = ifelse(length(ExpressionF()$Genes)>1,'Mean Expression', 'Expression'),
                 margin = list(l = 0,
                               r = 10,
                               b = 0,
                               t = 40,
                               pad = 0)) %>% 
          add_markers(x = ~reducedDim(sce,input$plotType)[,1], y=~reducedDim(sce,input$plotType)[,2],
                      z=~reducedDim(sce,input$plotType)[,3],
                      text= ~colData(sce)[,input$partitionType],
                      hoverinfo = 'text',
                      color = ~ExpressionF()$Exp,
                      name = ~colData(sce)[,input$partitionType],
                      size = I(point.size),span=I(0)) %>% 
          colorbar(title = "log(counts)",x=0,y=1) %>% 
          suppressWarnings(toWebGL())
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
                      text=  ~colData(sce)[,input$partitionType],
                      #name = ~colData(sce)[,input$partitionType],
                      hoverinfo = 'text') %>% 
          colorbar(title = "log(counts)") %>% 
          toWebGL()
      }
      
    })
    
    
    ### Heatmaps ----
    #Here I made everything in the same reactive object to manipulated all the reactive order at the same time
    Heatmap_Plot <- eventReactive(c(HeatmapF(),input$partitionType,input$mean_heat,input$cluster_row,input$cluster_column),{
      #req(input$scatter_heatmap == "heatmap")
      req(input$partitionType)
      req(!is.null(HeatmapF()))
      req(!is.null(OrderPartReact()))
      
      if(input$mean_heat) {
        dta <- apply(HeatmapF(),1,FUN =  function(x){
          tapply(x,colData(sce)[,input$partitionType],FUN=mean)
              }
        ) %>% t
        
        ht <-   HeatmapAnnotation(Type = levels(colData(sce)[,input$partitionType]),
                                  col=list(Type=OrderPartReact()$colPart),
                                  annotation_legend_param = list(Type = list(title = input$partitionType)
                                  ),
                                  annotation_label = c(input$partitionType),
                                  show_legend = c(Type =FALSE),
                                  show_annotation_name = T)
        
        h1 <- Heatmap(dta,
                      col = if(max(dta)==min(dta)) {viridis::viridis(1)} else {viridis::viridis(100)},
                      border =F,
                      name = "Gene expression",
                      cluster_rows = input$cluster_row,
                      cluster_columns = input$cluster_column,
                      row_names_side = "left",
                      column_names_rot = 45,
                      # column_title_rot = 45,
                      # column_title_gp = gpar(fontsize = 10),
                      column_title_side = "bottom",
                      # column_title = "Partition",
                      row_title = "Genes",
                      row_gap = unit(1, "mm"),
                      column_gap = unit(1, "mm"),
                      show_row_names = T,
                      # show_column_names = T,
                      top_annotation = ht,
                      # column_split = colData(sce)[,input$partitionType],
                      cluster_column_slices = F
                      # use_raster = TRUE,
                      # raster_by_magick = TRUE
        )
      } else{
        ht <-   HeatmapAnnotation(Type = colData(sce)[,input$partitionType],
                                  col=list(Type=OrderPartReact()$colPart),
                                  show_legend = c(Type =FALSE),
                                  annotation_label = c(input$partitionType),
                                  show_annotation_name = T)
        
        # colnames(HeatmapF()) <- NULL
        h1 <- Heatmap(as.matrix(HeatmapF()),
                      col = if(max(HeatmapF())==min(HeatmapF())) {viridis::viridis(1)} else {viridis::viridis(100)},
                      border =F,
                      name = "Gene expression",
                      cluster_rows = input$cluster_row,
                      cluster_columns = input$cluster_column,
                      row_names_side = "left",
                      # column_names_rot = 45,
                      column_title_rot = 45,
                      # column_title_gp = gpar(fontsize = 10),
                      column_title_side = "bottom",
                      # column_title = "Partition",
                      row_title = "Genes",
                      row_gap = unit(1, "mm"),
                      column_gap = unit(1, "mm"),
                      show_row_names = T,
                      show_column_names = F,
                      top_annotation = ht,
                      column_split = colData(sce)[,input$partitionType],
                      cluster_column_slices = F,
                      use_raster = TRUE,
                      raster_by_magick = TRUE
        )
      }
      h1
    })

    output$plot_heatmap <- renderPlot({
      #req(input$scatter_heatmap == "heatmap")
      req(!is.null(Heatmap_Plot()))
      Heatmap_Plot()
    })
    
    ####  Violin&SpikePlots ----
     output$Violin.Bar_Input <- renderUI({
       req(!is.null(ExpressionF()))
       radioGroupButtons(inputId = NS(id,"Cell_Exp"), label=NULL,choices = c("Violin", "SpikePlot"),
                         direction = "horizontal",justified = T,individual=T)
     })
     
     ViolinReact <- reactive({
       req(!is.null(ExpressionF()))
       data.frame(Y =ExpressionF()$Exp,
                  X=factor(colData(sce)[,input$partitionType]))
     })
     
     ViolinReact_Cell <- reactive({
       ViolinReact() %>% group_by(X) %>% summarise(n=n(),Ymax = (max(Y)+0.5))
     })
     
     ViolinPlot <-reactive({
       #req(input$scatter_heatmap == "scatter")
       req(!is.null(ExpressionF()))
       ggplot(ViolinReact()) + 
         geom_violin(aes(y = Y, 
                         x = X, 
                         fill = X),
                     scale="width") + 
         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
               legend.position = "none") + 
         xlab("Clusters") + 
         ylab("log(counts)") +
         ggtitle(ifelse(length(ExpressionF()$Genes)>1,"mean expression","expression")) +
         scale_fill_manual(values=OrderPartReact()$colPart) +
         geom_text(aes(label = n,x=X, y=Ymax), data = ViolinReact_Cell()) 
     })
     
     
     
     SpikePlot <-reactive({
       #req(input$scatter_heatmap == "scatter")
       req(!is.null(ExpressionF()))
       m <- barplot(ExpressionF()$Exp[OrderPartReact()$ordPart],
               col = OrderPartReact()$colPart[colData(sce)[,input$partitionType]][OrderPartReact()$ordPart],
               border = OrderPartReact()$colPart[colData(sce)[,input$partitionType]][OrderPartReact()$ordPart],
               ylab = "log(counts)", main = ifelse(length(ExpressionF()$Genes)>1,"mean expression","expression"), names.arg = F) 
		colleg <- legend_col(names(OrderPartReact()$colPart), max(m))
       legend(max(m)/2, -0.05, legend = names(OrderPartReact()$colPart), col = OrderPartReact()$colPart,
              pch=19, xpd=T, xjust = 0.5, cex = 0.9, ncol=colleg$ncol, text.width = colleg$colwidth)
       lines(x = m,
             tapply(ExpressionF()$Exp,
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
     
     
     #### Dotplots ----
     DotPlot  <- reactive({
        #req(input$scatter_heatmap == "dotplot")
        req(input$partitionType)
        if(input$GL_T){
          req(length(genes.GL()$genes) > 0)
          feature <- genes.GL()$genes
        } else { 
          req(!is.null(genes.L()))
          feature <- genes.L() 
        }
        
        g  <- scater::plotDots(object = sce,features = feature,group = input$partitionType,
                       scale = input$scale_dotplot,center = input$center_dotplot) + 
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=15)) +
          xlab(input$partitionType) + ylab("Genes")
        if(input$ord_dotplot){ #Define the order if the input is selected and the is not a heatmap of 1 row, received the order vector
          vtor <- if(length(feature)>1) {
            ord <- matrix(g$data$Average,byrow = F,nrow=length(unique(g$data$Feature))) %>% dist %>% hclust %>% .$order
            rev(feature[ord])
            } else {feature[1]}  
        } else {
          vtor <- rev(feature)
        }
        
         g$data$Feature <- factor(g$data$Feature,levels = vtor)
         g
         
         

     })
     
     output$plot_DotPlot <- renderPlotly({
       req(!is.null(DotPlot()))
       DotPlot() %>% ggplotly() %>% config(modeBarButtonsToRemove = c("select2d", "lasso2d", "hoverCompareCartesian"))
     })
     
     #### Stacked Violin ----
     stackVln <- reactive({
       #req(input$scatter_heatmap == "stackVln")
       req(input$partitionType)
       if(input$GL_T){
         req(length(genes.GL()$genes) > 0)
         feature <- genes.GL()$genes
       } else { 
         req(!is.null(genes.L()))
         feature <- genes.L() 
       }
       #Adaptated from https://github.com/ycl6/StackedVlnPlot
       df_plot <- assay(sce,"logcounts")[feature,,drop=F] %>% as.matrix %>% t %>% as.data.frame
       # Add cell ID and identity classes
       df_plot$Cell <- rownames(df_plot)
       df_plot$Idents <- colData(sce)[,input$partitionType]
       df_plot <- reshape2::melt(df_plot, id.vars = c("Cell","Idents"), measure.vars = feature,
                                 variable.name = "Feat", value.name = "Expr")
       
       if(input$ord_stackVln & length(feature) >1){ #Define the order if the input is selected and the is not a heatmap of 1 row, received the order vector
         dta <- apply(assay(sce,"logcounts")[feature,,drop=F],1,FUN =  function(x){
           tapply(x,colData(sce)[,input$partitionType],FUN=mean)
         }
         ) %>% t
         ord <- hclust(dist(dta))$order
         df_plot$Feat <- factor(df_plot$Feat, levels = levels(df_plot$Feat)[ord])
       } 
       
       
       g <- ggplot(df_plot, aes(x = factor(Idents),y =  Expr, fill = Idents)) +
         geom_violin(scale = "width", adjust = 1, trim = TRUE) +
         scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
           c(rep(x = "", times = max(length(x)-2, 0) ), x[length(x) - 1], "")) +
         facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
         cowplot::theme_cowplot(font_size = 12) +
         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
               legend.position = "none", panel.spacing = unit(0, "lines"),
               plot.title = element_text(hjust = 0.5),
               panel.background = element_rect(fill = NA, color = "black"),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold"),
               strip.text.y.left = element_text(angle = 0)) +
           xlab(input$partitionType) + ylab("Expression Level") + 
         scale_fill_manual(values=OrderPartReact()$colPart) 
       
       g
     })
     
     output$plot_stackVln <- renderPlot({
       req(!is.null(stackVln()))
       stackVln()
     })
     ### Downloads -----
     
     output$export_violin = downloadHandler(
       filename = function() {"Violin_Categories.pdf"},
       content = function(file) {
         pdf(file,
             width = input$pdf_width_violin,
             height = input$pdf_height_violin
         )
         ViolinPlot() %>% plot()
         dev.off()
       })
     
     output$export_SpikePlot = downloadHandler(
       filename = function() {"SpikePlot_Categories.pdf"},
       content = function(file) {
         pdf(file,
             width = input$pdf_width_SpikePlot,
             height = input$pdf_height_SpikePlot
         )
         SpikePlot() %>% print()
         dev.off()
       })
     
     output$export_heatmap = downloadHandler(
       filename = function() {"Heatmap_Categories.pdf"},
       content = function(file) {
         pdf(file,
             width = input$pdf_width_heatmap,
             height = input$pdf_height_heatmap
         )
         Heatmap_Plot() %>% plot()
         dev.off()
       }
     )
     
     output$export_dotplot = downloadHandler(
       filename = function() {"DotPlot_Categories.pdf"},
       content = function(file) {
         pdf(file,
             width = input$pdf_width_dotplot,
             height = input$pdf_height_dotplot
         )
         DotPlot() %>% plot()
         dev.off()
       }
     )
     
     output$export_StackedViolin = downloadHandler(
       filename = function() {"StackedViolin_Categories.pdf"},
       content = function(file) {
         pdf(file,
             width = input$pdf_width_StackedViolin,
             height = input$pdf_height_StackedViolin
         )
         stackVln() %>% plot()
         dev.off()
       }
     )
  }) 
}