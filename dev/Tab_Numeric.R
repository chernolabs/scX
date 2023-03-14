#####                              ######## 
#####        Partition Tab         ######## 
#####                              ######## 

##### Cluster UI Module ----
Clusters_UI <- function(id) {
  tagList(
    fluidRow(
      column(3,
             box(title = htmltools::span(icon("gears"), " Settings"),
                 width = NULL, status = "primary",solidHeader = T,collapsible = F,
                 conditionalPanel("input.scatter_heatmap == 'scatter'",ns=NS(id),
                   column(4, 
                          style='padding-left:12px; padding-right:3px;',
                          align="center",
                          pickerInput(NS(id,"partitionType1"),
                          "Partition #1",
                          choices = NULL)
                   ),
                   column(4,
                          style='padding-left:3px; padding-right:12px;',
                          align="center",
                          pickerInput(NS(id,"partitionType2"),
                          "Partition #2",
                          choices = NULL)
                   ),
                 ),
                 column(4,
                        style='padding-left:3px; padding-right:12px;',
                        align="center",
                        pickerInput(NS(id,"partitionColor"),
                                    "Partition Color",
                                    choices = NULL)
                 ),
                 conditionalPanel("input.scatter_heatmap == 'heatmap' || input.scatter_heatmap == 'dotplot' || input.scatter_heatmap == 'stackVln'", ns = NS(id),
                                  fluidRow(
                                    column(11,style='padding-left:0px; padding-right:2px;',
                                           selectizeInput(NS(id,"numeric_columns"),
                                                          label=NULL,
                                                          choices = NULL, 
                                                          options = list(maxItems = 10,
                                                                         maxOptions = 20,
                                                                         placeholder = 'Please select columns to plot'),
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
                                 plotlyOutput(NS(id,"plot_numeric"),height = "100vh") %>% withSpinner()
                             )
                    ),
                    tabPanel("Heatmap", value= "heatmap",
                             box(width = NULL, solidHeader = T, collapsible = F,
                                 footer = tagList(shiny::icon("cat"), "Nya"),
                                 plotOutput(NS(id,"plot_heatmap"),height = "100vh") %>% withSpinner()
                             )
                    ),
                    tabPanel("DotPlot", value= "dotplot",
                             box(width = NULL,solidHeader = T,collapsible = F,
                                 footer = tagList(shiny::icon("cat"), "Nya"),
                                 plotlyOutput(NS(id,"plot_DotPlot"),height = "100vh") %>% withSpinner()
                             )
                    ),
                    tabPanel("StackedViolin", value= "stackVln",
                             box(width = NULL,solidHeader = T,collapsible = F,
                                 footer = tagList(shiny::icon("cat"), "Nya"),
                                 plotOutput(NS(id,"plot_stackVln"),height = "100vh") %>% withSpinner()
                             )
                    )
             )
      )
    )
  )
}

##### Cluster Server Module ----
Clusters_Server <- function(id,sce) {
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
    
    output$plot_numeric <- renderPlotly({
      req(!is.null(df()))
      req(input$partitionType1)  
      req(input$partitionType2)
      req(input$partitionColor)  
      req(input$scatter_heatmap == "scatter")
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
       ggplotly(g) %>% toWebGL()
    })
    
    
    feature <- eventReactive(input$action,{
      req(length(input$numeric_columns) >0)
      input$numeric_columns
    })
    
    ### Heatmap ----
    Heatmap_DF <- eventReactive(c(feature(),input$norm_heat),{
      req(input$scatter_heatmap == "heatmap")
      req(length(feature()) > 0)
      
      # value <- ifelse(input$norm_heat,yes = "logcounts.norm",no = "logcounts") #To swtich between norm ot not normalize expression gene.
      exp_mtx <- df()[,feature(),drop=F] %>% as.matrix %>% t()
      
      if(input$norm_heat){
        exp_mtx <- apply(exp_mtx,1,function(x){x/max(x)}) %>% t()
      }
      
      exp_mtx
    })
    
    ### Heatmaps ----
    #Here I made everything in the same reactive object to manipulated all the reactive order at the same time
    Heatmap_Plot <- eventReactive(c(Heatmap_DF(),input$partitionColor,input$mean_heat,input$cluster_row,input$cluster_column),{
      req(input$scatter_heatmap == "heatmap")
      req(input$partitionColor)
      req(!is.null(Heatmap_DF()))
      
      if(input$mean_heat) {
        byPartition <- if(input$partitionColor != 'None'){df()[,input$partitionColor,drop=T]}else{rep('All',ncol(Heatmap_DF()))}
        dta <- apply(Heatmap_DF(),1,FUN =  function(x){
          tapply(x,byPartition,FUN=mean)
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
                      col = if(max(dta)==min(dta)) {viridis(1)} else {viridis(100)},
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
                      col = if(max(Heatmap_DF())==min(Heatmap_DF())) {viridis(1)} else {viridis(100)},
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
                      column_split = col_split,
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
    
    #### Dotplots ----
    output$plot_DotPlot <- renderPlotly({
      req(input$scatter_heatmap == "dotplot")
      req(input$partitionColor)
      req(length(feature()) > 0)
      
      byPartition <- if(input$partitionColor != 'None'){ input$partitionColor}else{NULL}
      g  <- plotDots(object = sce,features = feature(),group = ,
                     scale = input$scale_dotplot,center = input$center_dotplot) + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=15)) +
        xlab(input$partitionType) + ylab("Genes")
      if(input$ord_dotplot){ #Define the order if the input is selected and the is not a heatmap of 1 row, received the order vector
        vtor <- if(length(feature())>1) {
          ord <- matrix(g$data$Average,byrow = F,nrow=length(unique(g$data$Feature))) %>% dist %>% hclust %>% .$order
          rev(feature()[ord])
        } else {feature()[1]}  
      } else {
        vtor <- rev(feature())
      }
      
      g$data$Feature <- factor(g$data$Feature,levels = vtor)
      
      ggplotly(g) %>% config(modeBarButtonsToRemove = c("select2d", "lasso2d", "hoverCompareCartesian"))
      
      
    })
    
    
    ### Stck Violin
    output$plot_stackVln <- renderPlot({
      req(input$scatter_heatmap == "stackVln")
      req(length(feature()) > 0)
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
        theme_cowplot(font_size = 12) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              legend.position = "none", panel.spacing = unit(0, "lines"),
              plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = NA, color = "black"),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold"),
              strip.text.y.left = element_text(angle = 0)) +
         ylab("Expression Level")
      
        if(input$partitionColor == "None"){
          g <- g + scale_fill_manual(values="grey") + xlab("")
        } else{
          g <- g + scale_fill_manual(values=OrderPartReact()$colPart) + xlab(input$partitionColor)
        }
      g
    })
  
    
  })
}

ui <- fluidPage(
  Clusters_UI(id = 'lala')
)
server <- function(input, output, session) {
  Clusters_Server(sce = sce,id = 'lala')
}
shinyApp(ui = ui,server = server)
