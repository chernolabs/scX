#####                              ######## 
#####         Expression Tab       ######## 
#####                              ######## 

##### Expression UI Module ----
ExpressionUI <- function(id) {
  tagList(
    fluidRow(
      column(3,
        box(title = htmltools::span(icon("fa-light fa-gears"), " Settings"),
            width = NULL, status = "primary",solidHeader = T,collapsible = F,
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
          conditionalPanel("typeof output.plot !== 'undefined' || input.scatter_heatmap == 'heatmap' || input.scatter_heatmap == 'dotplot' || input.scatter_heatmap == 'stackVln'", ns = NS(id),
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
            prettySwitch(NS(id,"norm_heat"),
                         "Norm per gene",
                         value = F,
                         status = "primary",
                         fill = TRUE),
            prettySwitch(NS(id,"mean_heat"),
                         "Mean per group",
                         value = F,
                         status = "primary",
                         fill = TRUE),
            prettySwitch(NS(id,"cluster_row"),
                         "Cluster Row",
                         value = F,
                         status = "primary",
                         fill = TRUE),
            prettySwitch(NS(id,"cluster_column"),
                         "Cluster Column",
                         value = F,
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
          ),
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
                             icon = icon("fa-solid fa-play"))
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
                footer = tagList(shiny::icon("cat"), "Cat"),
              plotlyOutput(NS(id,"plot"),height = "100vh") %>% withLoader(type='html',loader = 'dnaspin')
            ),
            box(title = "Expression Plots",
                width = NULL, solidHeader = T,collapsible = T,
                footer = tagList(shiny::icon("cat"), "Cat"),
              uiOutput(NS(id,"Violin.Bar_Input")),
              plotOutput(NS(id,"Violin.Bar_Plot")) %>% withSpinner()
            )
          ),
          tabPanel("Heatmap", value= "heatmap",
            box(width = NULL, solidHeader = T, collapsible = F,
                footer = tagList(shiny::icon("cat"), "Cat"),
              plotOutput(NS(id,"plot_heatmap"),height = "100vh") %>% withLoader(type='html',loader = 'dnaspin')
            )
          ),
          tabPanel("DotPlot", value= "dotplot",
            box(width = NULL,solidHeader = T,collapsible = F,
                footer = tagList(shiny::icon("cat"), "Cat"),
              plotlyOutput(NS(id,"plot_DotPlot"),height = "100vh") %>% withLoader(type='html',loader = 'dnaspin')
            )
          ),
          tabPanel("StackedViolin", value= "stackVln",
            box(width = NULL,solidHeader = T,collapsible = F,
                footer = tagList(shiny::icon("cat"), "Cat"),
              plotOutput(NS(id,"plot_stackVln"),height = "100vh") %>% withLoader(type='html',loader = 'dnaspin')
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
    
    updatePickerInput(session, 'partitionType', 
                         choices = names(colData(sce))[sapply(colData(sce), is.factor)])
    
    updateSelectizeInput(session, 'gen_exp', choices = rownames(sce), server = TRUE)
    
    dimVector <- reactive({
        sapply(reducedDims(x = sce),FUN = ncol) 
    })
    
    observeEvent(c(dimVector()),{
      req(input$scatter_heatmap == "scatter")
      updatePickerInput(session,inputId = "DimType", choices = (c("3","2")[c(3,2) %in% dimVector()]))
    })
    
    observeEvent(c(input$DimType,dimVector()), {
      req(input$scatter_heatmap == "scatter")
      req(!is.null(dimVector()))
      req(input$DimType)
	  updatePickerInput(session,inputId = "plotType", choices = rev(names(which(dimVector() == as.numeric(input$DimType)  | dimVector() > 3))))
      # if(input$DimType == "3"){
        # updatePickerInput(session,inputId = "plotType", choices = rev(names(which(dimVector() == as.numeric(input$DimType)  | dimVector() > 3))))
      # } else if(input$DimType == "2") { 
        # updatePickerInput(session,inputId = "plotType", choices = rev(names(which(dimVector() == as.numeric(input$DimType)  | dimVector() > 3))))
      # }
    })
    
    ### Gen selected & data preparation ----
        #### Selected Gene -----
    genes.L <- eventReactive(input$action,{
      req(input$GL_T == F)
      req(input$gen_exp)
      input$gen_exp
    })
    
    ExpressionL <- eventReactive(genes.L(),{
      req(input$scatter_heatmap == "scatter")
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
      req(input$scatter_heatmap == "heatmap")
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
      req(input$scatter_heatmap == "scatter")
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
      req(input$scatter_heatmap == "heatmap")
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
      req(input$scatter_heatmap == "heatmap")
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
    ClusterPlot <- eventReactive(c(input$plotType,input$partitionType),{
      req(input$scatter_heatmap == "scatter")
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
                      size = I(point.size),span=I(0),text=~colData(sce)[,input$partitionType],hoverinfo='text') %>% 
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
                      size = I(point.size),span=I(0),text=~colData(sce)[,input$partitionType],hoverinfo='text') %>% 
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
                      text=  ~colData(sce)[,input$partitionType],
                      #name = ~colData(sce)[,input$partitionType],
                      hoverinfo = 'text') %>% 
          colorbar(title = "log(counts)") %>% 
          toWebGL()
      }
      
    })
    
    output$plot <- renderPlotly({
      if(input$button == T){
        ClusterPlot()
      } else{
        ExpressionPlot()
      }
    })
    
    
    
    ### Heatmaps ----
    #Here I made everything in the same reactive object to manipulated all the reactive order at the same time
    Heatmap_Plot <- eventReactive(c(HeatmapF(),input$partitionType,input$mean_heat,input$cluster_row,input$cluster_column),{
      req(input$scatter_heatmap == "heatmap")
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
                      col = viridis(100),
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
                      col = viridis(100),
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
      req(input$scatter_heatmap == "heatmap")
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
       req(!is.null(ExpressionF()) & input$Cell_Exp == "Violin")
       data.frame(Y =ExpressionF()$Exp,
                  X=factor(colData(sce)[,input$partitionType]))
     })
     
     ViolinReact_Cell <- reactive({
       ViolinReact() %>% group_by(X) %>% summarise(n=n(),Ymax = (max(Y)+0.5))
     })
     
     ViolinPlot <-reactive({
       req(input$scatter_heatmap == "scatter")
       req(!is.null(ExpressionF()) & input$Cell_Exp == "Violin")
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
       req(input$scatter_heatmap == "scatter")
       req(!is.null(ExpressionF()) & input$Cell_Exp == "SpikePlot")
       barplot(ExpressionF()$Exp[OrderPartReact()$ordPart],
               col = OrderPartReact()$colPart[colData(sce)[,input$partitionType]][OrderPartReact()$ordPart],
               border = OrderPartReact()$colPart[colData(sce)[,input$partitionType]][OrderPartReact()$ordPart],
               ylab = "log(counts)", main = ifelse(length(ExpressionF()$Genes)>1,"mean expression","expression"), names.arg = F) 
       legend("bottom", legend = names(OrderPartReact()$colPart), col = OrderPartReact()$colPart,
              pch=19, ncol=6, xpd=T, inset=c(0,-0.25))
       abline(h=mean(ExpressionF()$Exp[ExpressionF()$Exp>0]),lty=2,col="grey")
     })
     
     output$Violin.Bar_Plot <- renderPlot({
       req(!is.null(ExpressionF()))
       req(input$Cell_Exp)
       if(input$Cell_Exp == "Violin"){
         ViolinPlot()
       } else { 
         SpikePlot()
       }
     })
     
     
     #### Dotplots ----
     output$plot_DotPlot <- renderPlotly({
        req(input$scatter_heatmap == "dotplot")
        req(input$partitionType)
        if(input$GL_T){
          req(length(genes.GL()$genes) > 0)
          feature <- genes.GL()$genes
        } else { 
          req(!is.null(genes.L()))
          feature <- genes.L() 
        }
        
        g  <- plotDots(object = sce,features = feature,group = input$partitionType,
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
         
         ggplotly(g)
         

     })
     
     #### Stacked Violin ----
     output$plot_stackVln <- renderPlot({
       req(input$scatter_heatmap == "stackVln")
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
           c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
         facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
         theme_cowplot(font_size = 12) +
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
  }) 
}