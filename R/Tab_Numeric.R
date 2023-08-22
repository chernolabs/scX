#####                              ######## 
#####        Partition Tab         ######## 
#####                              ######## 

##### Cluster UI Module ----
Fields_UI <- function(id) {
  tagList(
    fluidRow(
      column(3,
             box(title = htmltools::span(icon("gears"), " Settings"),
                 width = NULL, status = "primary",solidHeader = T,collapsible = F,
                 conditionalPanel("input.scatter_heatmap == 'scatter'",ns=NS(id),
                   fluidRow(
                     column(6, 
                            style='padding-left:12px; padding-right:3px;',
                            align="center",
                            pickerInput(NS(id,"partitionType1"),
                            "Field #1",
                            choices = NULL)
                     ),
                     column(6,
                            style='padding-left:3px; padding-right:12px;',
                            align="center",
                            pickerInput(NS(id,"partitionType2"),
                            "Field #2",
                            choices = NULL)
                     )
                   ),
                 ),
                 conditionalPanel("input.scatter_heatmap == 'heatmap' || input.scatter_heatmap == 'dotplot' || input.scatter_heatmap == 'stackVln' || input.scatter_heatmap == 'matrix' ", ns = NS(id),
                                  fluidRow(
                                    column(10,style='padding-left:12px; padding-right:2px;',align="center",
                                           selectizeInput(NS(id,"numeric_columns"),
                                                          label=NULL,
                                                          choices = NULL, 
                                                          options = list(maxItems = 10,
                                                                         maxOptions = 20,
                                                                         placeholder = 'Please select columns to plot'),
                                                          width = NULL,
                                                          multiple=T)
                                    ),
                                    column(2,style='padding-left:2px; padding-right:12px; padding-top:4px',
                                           align="center",
                                           actionBttn((NS(id,"action")),
                                                      label = NULL,
                                                      style = "unite",
                                                      color = "primary",
                                                      size = "xs",
                                                      icon = icon("play"))
                                    )
                                  )
                 ),
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
                                               "Norm per field",
                                               value = F,
                                               status = "primary",
                                               fill = TRUE),
                                  prettySwitch(NS(id,"mean_heat"),
                                               "Mean per group",
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
                 conditionalPanel("input.scatter_heatmap == 'matrix'",ns=NS(id),
                                  prettySwitch(NS(id,"matrix_cluster_row"),
                                               "Cluster Row",
                                               value = F,
                                               status = "primary",
                                               fill = TRUE),
                                  prettySwitch(NS(id,"matrix_cluster_column"),
                                               "Cluster Column",
                                               value = F,
                                               status = "primary",
                                               fill = TRUE),
                                  prettySwitch(NS(id,"showFreq"),
                                               "Hide Labels",
                                               value = F,
                                               status = "primary",
                                               fill = TRUE)
                                  
                 ),
				 hr(style = "border-top: 1px solid #0073b7;"),
				 column(12,
                        style='padding-left:20px; padding-right:20px;',
                        align="center",
                        pickerInput(NS(id,"partitionColor"),
                                    "Partition Color",
                                    choices = NULL)
                 )
             )
      ),
      column(9,
             tabBox(id = NS(id,"scatter_heatmap"),
                    selected = "scatter",
                    width = NULL,
                    tabPanel("Distribution",value = "scatter",
                             box(title = "Distribution Plot",
                                 width = NULL, solidHeader = T, collapsible = T,
                                 footer = tagList(shiny::icon("cat"), "Nya"),
                                 dropdownButton(
									fluidRow(
										column(4, style='padding-left:12px; padding-right:3px;', numericInput(NS(id,"pdf_width_scatter"),"Width",value = 7)),
										column(4, style='padding-left:3px; padding-right:3px;', numericInput(NS(id,"pdf_height_scatter"),"Height",value = 7)),
										column(4, style='padding-left:3px; padding-right:12px; padding:16px', downloadButton(NS(id,'export_scatter')))
									),
                                   circle = FALSE,
                                   status = "primary",
                                   icon = icon("fa-thin fa-download"),
                                   width = "300px",
                                   size= "sm",
                                   up = F,
                                   tooltip = tooltipOptions(title = "Press to Download")
                                 ),
                                 plotlyOutput(NS(id,"plot_numeric"),height = "80vh") %>% withSpinner()
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
                                 plotOutput(NS(id,"plot_heatmap"),height = "80vh") %>% withSpinner()
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
                                   icon = icon("fa-thin fa-download"),
                                   width = "300px",
                                   size= "sm",
                                   up = F,
                                   tooltip = tooltipOptions(title = "Press to Download")
                                 ),
                                 plotlyOutput(NS(id,"plot_DotPlot"),height = "80vh") %>% withSpinner()
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
                                   icon = icon("fa-thin fa-download"),
                                   width = "300px",
                                   size= "sm",
                                   up = F,
                                   tooltip = tooltipOptions(title = "Press to Download")
                                 ),
                                 plotOutput(NS(id,"plot_stackVln"),height = "80vh") %>% withSpinner()
                             )
                    ),
                    tabPanel("Matrix", value= "matrix",
                                 uiOutput(NS(id,"CorheatMapOutput"))
                    )
             )
      )
    )
  )
}

##### Cluster Server Module ----
Fields_Server <- function(id,sce) {
  moduleServer(id, function(input,output,session) {
    
    ### Observe Events ----
    updatePickerInput(session, 'partitionColor', 
                      choices = c('None',names(colData(sce))[sapply(colData(sce), is.factor)]))
    
    updatePickerInput(session,inputId = "partitionType1", 
                      choices = names(colData(sce))[sapply(colData(sce), is.numeric)])
    
    updatePickerInput(session,inputId = "partitionType1", 
                        choices = names(colData(sce))[sapply(colData(sce), is.numeric)])
    
    observeEvent(input$partitionType1,ignoreInit = T,{
      req(input$partitionType1)
      prt <- names(colData(sce))[sapply(colData(sce), is.numeric)]
      
      updatePickerInput(session,inputId = "partitionType2", 
                        choices = c("None",prt[-(match(input$partitionType1,prt))])
      )
    })
    
    updateSelectizeInput(session, 'numeric_columns', choices =names(colData(sce))[sapply(colData(sce), is.numeric)], server = TRUE)
    
    # observeEvent(input$partitionType2, {
    #   if(input$partitionType2 == "None"){
    #     hideTab(inputId = "scatter_heatmap", target = "matrix")  
    #   } else {
    #     showTab(inputId = "scatter_heatmap", target = "matrix")
    #   }
    # })
    
    #### Plots ----
    
    OrderPartReact <- eventReactive(input$partitionColor,{
      req(input$partitionColor != 'None') 
      Col.and.Order(partition = input$partitionColor, sce=sce)
    })
    
    df <- reactive({
        colData(sce) %>% as.data.frame()
    })
    
    ### scatters ----
    
    PlotNumeric <- reactive({
      req(!is.null(df()))
      req(input$partitionType1)  
      req(input$partitionType2)
      req(input$partitionColor)  
      #req(input$scatter_heatmap == "scatter")
      if(input$partitionType2 == "None"){
        if(input$partitionColor == "None"){
          g <- df() %>% ggplot(aes(y=.data[[input$partitionType1]],x='All',fill='All')) + geom_boxplot() +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none") + 
            ylab(input$partitionType1) + scale_fill_manual(values="grey")
        } else{
          g <- df() %>% ggplot(aes(y=.data[[input$partitionType1]],x=.data[[input$partitionColor]],fill=.data[[input$partitionColor]])) + geom_boxplot() +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none") + 
            ylab(input$partitionType1) +
            xlab(input$partitionColor) + scale_fill_manual(values=OrderPartReact()$colPart)
        }
      }
      else{
        if(input$partitionColor == "None"){
          g <- df() %>% ggplot(aes(x=.data[[input$partitionType1]],y=.data[[input$partitionType2]])) + geom_point(alpha=0.5,colour="black") +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none") + 
             xlab(input$partitionType1) + ylab(input$partitionType2)
        } else{
          g <- df() %>% ggplot(aes(x=.data[[input$partitionType1]],y=.data[[input$partitionType2]],col=.data[[input$partitionColor]])) + geom_point(alpha=0.5) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
            xlab(input$partitionType1) + ylab(input$partitionType2) + scale_color_manual(values=OrderPartReact()$colPart)
        }
      }
       g
    })
    
    output$plot_numeric <- renderPlotly({
      req(!is.null(PlotNumeric()))
      PlotNumeric() %>% ggplotly() %>% toWebGL()
      })
    
    feature <- eventReactive(input$action,{
      req(length(input$numeric_columns) >0)
      input$numeric_columns
    })
    
    ### Heatmap ----
    Heatmap_DF <- eventReactive(c(feature(),input$norm_heat),{
      #req(input$scatter_heatmap == "heatmap")
      req(length(feature()) > 0)
      
      # value <- ifelse(input$norm_heat,yes = "logcounts.norm",no = "logcounts") #To swtich between norm ot not normalize expression gene.
      exp_mtx <- df()[,feature(),drop=F] %>% as.matrix %>% t()
      
      if(input$norm_heat){
        exp_mtx <- apply(exp_mtx,1,function(x){x/max(x,na.rm = T)}) %>% t()
      }
      
      exp_mtx
    })
    
    ### Heatmaps ----
    #Here I made everything in the same reactive object to manipulated all the reactive order at the same time
    Heatmap_Plot <- eventReactive(c(Heatmap_DF(),input$partitionColor,input$mean_heat,input$cluster_row,input$cluster_column),{
      #req(input$scatter_heatmap == "heatmap")
      req(input$partitionColor)
      req(!is.null(Heatmap_DF()))
      
      if(input$mean_heat) {
        byPartition <- if(input$partitionColor != 'None'){df()[,input$partitionColor,drop=T]}else{rep('All',ncol(Heatmap_DF()))}
        dta <- apply(Heatmap_DF(),1,FUN =  function(x){
          tapply(x,byPartition,FUN=function(y){mean(y,na.rm = T)})
        }
        )
        
        if(input$partitionColor != 'None'){
        dta <- dta %>% t
        ht <-   HeatmapAnnotation(Type = levels(df()[,input$partitionColor]),
                                  col=list(Type=OrderPartReact()$colPart),
                                  annotation_legend_param = list(Type = list(title = input$partitionColor)
                                  ),
                                  annotation_label = c(input$partitionColor),
                                  show_legend = c(Type =FALSE),
                                  show_annotation_name = T)
        } else{
          ht <- NULL
        }
        
        h1 <- Heatmap(dta,
                      col = if(max(dta,na.rm = T)==min(dta,na.rm = T)) {viridis::viridis(1)} else {viridis::viridis(100)},
                      border =F,
                      na_col = "grey",
                      name = "Value",
                      cluster_rows = input$cluster_row,
                      cluster_columns = input$cluster_column,
                      row_names_side = "left",
                      column_names_rot = 45,
                      # column_title_rot = 45,
                      # column_title_gp = gpar(fontsize = 10),
                      column_title_side = "bottom",
                      # column_title = "Partition",
                      row_title = "Fields",
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
        if(input$partitionColor != 'None'){
        ht <-   HeatmapAnnotation(Type = df()[,input$partitionColor],
                                  col=list(Type=OrderPartReact()$colPart),
                                  show_legend = c(Type =FALSE),
                                  annotation_label = c(input$partitionColor),
                                  show_annotation_name = T)
        col_split <- df()[,input$partitionColor]
        } else{
          ht <- NULL
          col_split <- NULL
        }
        # colnames(HeatmapF()) <- NULL
        h1 <- Heatmap(as.matrix(Heatmap_DF()),
                      col = if(max(Heatmap_DF(),na.rm = T)==min(Heatmap_DF(),na.rm = T)) {viridis::viridis(1)} else {viridis::viridis(100)},
                      border =F,
                      na_col = "grey",
                      name = "Value",
                      cluster_rows = input$cluster_row,
                      cluster_columns = input$cluster_column,
                      row_names_side = "left",
                      # column_names_rot = 45,
                      column_title_rot = 45,
                      # column_title_gp = gpar(fontsize = 10),
                      column_title_side = "bottom",
                      # column_title = "Partition",
                      row_title = "Fields",
                      row_gap = unit(1, "mm"),
                      column_gap = unit(1, "mm"),
                      show_row_names = T,
                      show_column_names = F,
                      top_annotation = ht,
                      column_split = col_split,
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
    
    #### Dotplots ----
    DotPlot <- reactive({
      req(length(feature()) > 0)
      #req(input$scatter_heatmap == "dotplot")
      req(input$partitionColor)
      
      byPartition <- if(input$partitionColor != 'None'){ input$partitionColor}else{NULL}
      g  <- plotDots_fields(df = df(),
                            features = feature(),
                            group = byPartition,
                     scale = input$scale_dotplot,center = input$center_dotplot) + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=15)) +
        xlab(input$partitionType) + ylab("Fields")
      if(input$ord_dotplot){ #Define the order if the input is selected and the is not a heatmap of 1 row, received the order vector
        vtor <- if(length(feature())>1) {
          ord <- matrix(g$data$Average,byrow = F,nrow=length(unique(g$data$Feature))) %>% dist %>% hclust %>% .$order
          rev(feature()[ord])
        } else {feature()[1]}  
      } else {
        vtor <- rev(feature())
      }
      
      g$data$Feature <- factor(g$data$Feature,levels = vtor)
      
      g
    })
    
    output$plot_DotPlot <- renderPlotly({
      req(!is.null(DotPlot()))
      DotPlot() %>% ggplotly() %>% config(modeBarButtonsToRemove = c("select2d", "lasso2d", "hoverCompareCartesian"))
    })
    
    ### Stck Violin
    stackVln <- reactive({
      req(length(feature()) > 0)
      #req(input$scatter_heatmap == "stackVln")
      #Adaptated from https://github.com/ycl6/StackedVlnPlot
      df_plot <- df()[,feature(),drop=F] %>% as.data.frame
      # Add cell ID and identity classes
      df_plot$Cell <- rownames(df_plot)
      
      if(input$partitionColor == "None") {
        df_plot$Idents <- "All"
        df_plot$Idents <- as.factor(df_plot$Idents)  
      } else {
        df_plot$Idents <- colData(sce)[,input$partitionColor]
      }
      
      df_plot <- reshape2::melt(df_plot, id.vars = c("Cell","Idents"), measure.vars = feature(),
                                variable.name = "Feat", value.name = "Expr")
      
      if(input$ord_stackVln & length(feature()) >1){ #Define the order if the input is selected and the is not a heatmap of 1 row, received the order vector
        dta <- apply(t(as.matrix(colData(sce)[,feature(),drop=F])),1,FUN =  function(x){
          tapply(X = x,
                 INDEX = if(partitionColor == "None"){
                          rep("All",ncol(sce))
                         } else{
                           colData(sce)[,partitionColor]
                         },
                 FUN=mean)
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
         ylab("Value")
      
        if(input$partitionColor == "None"){
          g <- g + scale_fill_manual(values="grey") + xlab("")
        } else{
          g <- g + scale_fill_manual(values=OrderPartReact()$colPart) + xlab(input$partitionColor)
        }
      g
    })
    
    output$plot_stackVln <- renderPlot({
      req(!is.null(stackVln()))
      stackVln()
    })
    ### Matrix ----
    
    Matrix_DF <- eventReactive(c(input$partitionColor,feature()),{
      #req(input$scatter_heatmap == "matrix")
      req(length(feature()) > 1)
      byPartition <- if(input$partitionColor != 'None'){input$partitionColor}else{NULL}
      df_corr <- df() %>% group_by(across(all_of(byPartition))) %>% summarise(Correlation = cor(cbind(across(all_of(feature())))),use = "complete.obs",method = "spearman")
      df_corr
    })
    
    PlotMatrix <- reactive({
      #req(input$scatter_heatmap == "matrix")
      req(!is.null(Matrix_DF()))
      if(input$partitionColor != 'None'){
        # rownames(df_corr$Correlation)  <- paste0(mtx[,input$partitionColor,drop=T],"_",rownames(mtx$Correlation))
        ht <-   HeatmapAnnotation(Type = Matrix_DF()[,input$partitionColor,drop=T],
                                  col=list(Type=OrderPartReact()$colPart),
                                  annotation_legend_param = list(Type = list(title = input$partitionColor)
                                  ),
                                  annotation_label = c(input$partitionColor),
                                  show_legend = c(Type =FALSE),
                                  show_annotation_name = T)
        col_split <- Matrix_DF()[,input$partitionColor,drop=T]
        mtx <- Matrix_DF()$Correlation %>% t()
        rownames(mtx) <- feature()
      } else{
        mtx <- Matrix_DF()$Correlation
        ht <- NULL
        col_split <- NULL
      }
      Heatmap(mtx,
              col = if(max(t(mtx),na.rm = T)==min(t(mtx),na.rm = T)) {viridis::viridis(1)} else {viridis::viridis(100)},
              border =F,
              name = "Correlation",
              na_col = "grey",
              cluster_rows = input$matrix_cluster_row,
              cluster_columns = input$matrix_cluster_column,
              row_names_side = "left",
              column_title_rot = 45,
              column_title_side = "bottom",
              row_title = "Fields",
              row_gap = unit(1, "mm"),
              column_gap = unit(1, "mm"),
              show_row_names = T,
              show_column_names = T,
              top_annotation = ht,
              column_split = col_split,
              cluster_column_slices = F,
              cell_fun = if(!(input$showFreq)){set_val(tab = mtx)} else{NULL}
      )
    })
    
    output$plot_Matrix <- renderPlot({
      req(!is.null(PlotMatrix()))
      PlotMatrix()
    })
    
    output$CorheatMapOutput <- renderUI({
      #req(input$scatter_heatmap == "matrix")
      req(!is.null(Matrix_DF()))
      #If there are only 1 and NA values, it doens't show any plot.
      if(all(unique(as.vector(Matrix_DF()$Correlation))%in% c(NA,1))){
        HTML('<p style="text-align: center;"><strong>No correlation could be found with this conditions.</strong></p>')  
      } else{
        tagList(
          box(width = NULL, solidHeader = T, collapsible = F,
              footer = tagList(shiny::icon("cat"), "Nya"),
              dropdownButton(
				fluidRow(
					column(4, style='padding-left:12px; padding-right:3px;', numericInput(NS(id,"pdf_width_matrix"),"Width",value = 7)),
					column(4, style='padding-left:3px; padding-right:3px;', numericInput(NS(id,"pdf_height_matrix"),"Height",value = 7)),
					column(4, style='padding-left:3px; padding-right:12px; padding:16px', downloadButton(NS(id,'export_matrix')))
				),
                circle = FALSE,
                status = "primary",
                icon = icon("fa-thin fa-download"),
                width = "300px",
                size= "sm",
                up = F,
                tooltip = tooltipOptions(title = "Press to Download")
              ),
        plotOutput(NS(id,"plot_Matrix"),height = "80vh") %>% withSpinner()
          )
        )
      }
    })
    
    ### Downloads -----
    
    output$export_scatter = downloadHandler(
      filename = function() {"Distribution_Numerics.pdf"},
      content = function(file) {
        pdf(file,
            width = input$pdf_width_scatter,
            height = input$pdf_height_scatter
        )
        PlotNumeric() %>% plot()
        dev.off()
      })
    
    output$export_heatmap = downloadHandler(
      filename = function() {"Heatmap_Numerics.pdf"},
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
      filename = function() {"DotPlot_Numerics.pdf"},
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
      filename = function() {"StackedViolin_Numerics.pdf"},
      content = function(file) {
        pdf(file,
            width = input$pdf_width_StackedViolin,
            height = input$pdf_height_StackedViolin
        )
        stackVln() %>% plot()
        dev.off()
      }
    )
    
    output$export_matrix = downloadHandler(
      filename = function() {"CorrMatrix_Numerics.pdf"},
      content = function(file) {
        pdf(file,
            width = input$pdf_width_matrix,
            height = input$pdf_height_matrix
        )
        PlotMatrix() %>% plot()
        dev.off()
      }
    )
    
  })
}
# 
# ui <- fluidPage(
#   Clusters_UI(id = 'lala')
# )
# server <- function(input, output, session) {
#   Clusters_Server(sce = sce,id = 'lala')
# }
# shinyApp(ui = ui,server = server)

