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
    
    
    output$plot_stackVln <- renderPlot({
      req(input$scatter_heatmap == "stackVln")
      req(length(feature()) > 0)
      #Adaptated from https://github.com/ycl6/StackedVlnPlot
      df_plot <- colData(sce)[,feature(),drop=F] %>% as.data.frame
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
