#####                              ######## 
#####        Partition Tab         ######## 
#####                              ######## 

##### Cluster UI Module ----
Clusters_UI <- function(id) {
  tagList(
    fluidRow(
      column(3,
        box(title = htmltools::span(icon("gears"), " Settings"),
            width = NULL, status = "primary",solidHeader = T,collapsible = T,
          fluidRow(
            column(6, 
                   style='padding-left:12px; padding-right:3px;',
                   align="center",
              pickerInput(NS(id,"partitionType1"),
                          "Partition #1",
                          choices = NULL)),
            column(6,
                   style='padding-left:3px; padding-right:12px;',
                   align="center",
              pickerInput(NS(id,"partitionType2"),
                          "Partition #2",
                          choices = NULL)
            )
          ),
		  fluidRow(
			column(6, 
				conditionalPanel("(input.barplot_matrix == 'barplot' && (input.wrap || input.partitionType2 == 'None')) || input.barplot_matrix == 'matrix'", ns=NS(id),  
					prettySwitch(NS(id,"showFreq"),
                       label="Hide Labels",
                       value = F,
                       status = "primary",
                       fill = TRUE
					)
				),
				conditionalPanel("input.barplot_matrix == 'barplot' && (input.wrap || input.partitionType2 == 'None')",ns=NS(id),
					prettySwitch(NS(id,"freq"),
                           label= "Proportions",
                           value = F,
                           status = "primary",
                           fill = TRUE
					)
				)
			),
			column(6, 
				conditionalPanel("input.barplot_matrix == 'barplot' && input.partitionType2 != 'None'",ns=NS(id),
					prettySwitch(NS(id,"wrap"),
                           label="Wrap",
                           value = F,
                           status = "primary",
                           fill = TRUE
					)
				)
			)
		  ),
          conditionalPanel("input.barplot_matrix == 'matrix' && input.partitionType2 != 'None'",ns=NS(id),
				htmlOutput(NS(id,"randValue"), style = "text-align:center;"),
				bsTooltip(NS(id, "randValue"), title = "A measure of similarity between two partitions", placement = "bottom", trigger = "hover", options = list(container = "body")),
             radioGroupButtons(inputId = NS(id,"Metric"), 
                               label="",
                               choices = c("Count", "Jaccard"),
                               direction = "horizontal",
                               justified = T,
                               individual=F)
          )
        )
      ),
      column(9,
        tabBox(id = NS(id,"barplot_matrix"),
               selected = "barplot",
               width = NULL,
          tabPanel("BarPlot",value = "barplot",
                   dropdownButton(
					fluidRow(
						column(4, style='padding-left:12px; padding-right:3px;', numericInput(NS(id,"pdf_width_barplot"),"Width",value = 7)),
						column(4, style='padding-left:3px; padding-right:3px;', numericInput(NS(id,"pdf_height_barplot"),"Height",value = 7)),
						column(4, style='padding-left:3px; padding-right:12px; padding:16px', downloadButton(NS(id,'export_barplot')))
					),
                     circle = FALSE,
                     status = "primary",
                     icon = icon("download"),
                     width = "300px",
                     size= "sm",
                     up = F,
                     tooltip = tooltipOptions(title = "Download")
                   ),
            plotOutput(NS(id,"plot_cluster"),
                     height = "80vh") %>% withSpinner()
          ),
          tabPanel("Matrix",value = "matrix",
                   dropdownButton(
					fluidRow(
						column(4, style='padding-left:12px; padding-right:3px;', numericInput(NS(id,"pdf_width_matrix"),"Width",value = 7)),
						column(4, style='padding-left:3px; padding-right:3px;', numericInput(NS(id,"pdf_height_matrix"),"Height",value = 7)),
						column(4, style='padding-left:3px; padding-right:12px; padding:16px', downloadButton(NS(id,'export_matrix')))
					),
                     circle = FALSE,
                     status = "primary",
                     icon = icon("download"),
                     width = "300px",
                     size= "sm",
                     up = F,
                     tooltip = tooltipOptions(title = "Download")
                   ),
            plotOutput(NS(id,"plot_matrix"),
                       height = "80vh") %>% withSpinner()
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
    
    updatePickerInput(session,inputId = "partitionType1", 
                        choices = names(colData(sce))[sapply(colData(sce), is.factor)])
    
    observeEvent(input$partitionType1,ignoreInit = T,{
      req(input$partitionType1)
      prt <- names(colData(sce))[sapply(colData(sce), is.factor)]
      
      updatePickerInput(session,inputId = "partitionType2", 
                        choices = c("None",prt[-(match(input$partitionType1,prt))])
      )
    })
    
    observeEvent(input$partitionType2, {
      if(input$partitionType2 == "None"){
        hideTab(inputId = "barplot_matrix", target = "matrix")  
      } else {
        showTab(inputId = "barplot_matrix", target = "matrix")
      }
    })
    
    #### Plots ----
    OrderPartReact <- eventReactive(c(input$partitionType1,
                                      input$partitionType2),{
      req(input$partitionType1)
      req(input$partitionType2)
      if(input$partitionType2 == "None"){
        Col.and.Order(partition = input$partitionType1, sce=sce)
      } else {
        Col.and.Order(partition = input$partitionType2, sce=sce)
      }
    })
    
    BarplotData <- reactive({
      if(input$partitionType2 == "None"){
      df <- data.frame(Partition1 = colData(sce)[,input$partitionType1]) %>%
          group_by(Partition1) %>% summarise(Val =n()) %>% 
          mutate(Freq = Val / sum(Val)) %>% as.data.frame
      df
      } else {
      df <- data.frame(Partition1 = colData(sce)[,input$partitionType1],
                   Partition2 = colData(sce)[,input$partitionType2]) %>%
            group_by(Partition1,Partition2) %>% 
            summarise(Val = n()) 
      
      }
      df$Partition1 <- factor(df$Partition1,levels=rev(levels(df$Partition1)))
      df
    })
    
    ### BarPlots ----
    
    Barplot_cluster <- reactive({
      req(input$partitionType1)  
      req(input$partitionType2)  
      # req(input$barplot_matrix == "barplot")
      if(input$partitionType2 == "None"){
        ### Simple BarPlot ----
              ### Proportion ----
        if(input$freq){
          g  <- ggplot(BarplotData(),
                       aes(x= Partition1,
                           fill= Partition1,
                           y= Freq)) +
                       geom_bar(stat = "identity") + 
                       coord_flip() +
                       scale_fill_manual(values=OrderPartReact()$colPart) +
                       theme(legend.position = "none") +
                       ylab("Relative Frequencies") + 
                       xlab(input$partitionType1)
          if(!(input$showFreq)){
            g <- g + geom_text( aes(label=round(Freq,2)), position = position_stack(vjust = 0.85))
          }
        } else {
              #### Counts ----
          g  <-  ggplot(BarplotData(),
                        aes(x=Partition1,
                            fill=Partition1,
                            y= Val)) +
           geom_bar(stat = "identity") + 
                       coord_flip() +
                       scale_fill_manual(values=OrderPartReact()$colPart) +
                       theme(legend.position = "none") +
                       ylab("Absolute Frequencies") + 
                       xlab(input$partitionType1)
          if(!(input$showFreq)){
            g <- g + geom_text( aes(label=round(Val,2)), position = position_stack(vjust = 0.85))
          }
        }
        ### 2 partitons BarPlot ----
      } else {
          if(input$wrap) {
            if(input$freq){
              df <- BarplotData() %>% group_by(Partition2) %>% 
                              mutate(Freq = Val / sum(Val))
              
              g  <- ggplot(df,aes(x = Partition1,
                                  fill=Partition2,
                                  y = Freq)) +
                geom_bar(stat = "identity") +
                facet_grid(.~Partition2) +
                coord_flip() + 
                xlab(input$partitionType1) + 
                ylab("Proportion") +
                scale_fill_manual(values=OrderPartReact()$colPart) +
                theme(axis.text.x = element_text(angle = 45,
                                                 vjust = 1,
                                                 hjust=1),
                      legend.position = "none")
              if(!(input$showFreq)){
                g <- g + geom_text( aes(label=round(Freq,2),y = (max(Freq)/2)))
              }
            } else {
              g  <- ggplot(BarplotData(),
                     aes(x = Partition1,
                         fill=Partition2,
                         y = Val)) +
                  geom_bar(stat = "identity") +
                  facet_grid(.~Partition2) +
                  coord_flip() + 
                  scale_fill_manual(values=OrderPartReact()$colPart)  +
                  theme(axis.text.x = element_text(angle = 45,
                                                   vjust = 1,
                                                   hjust=1),
                        legend.position = "none")
              if(!(input$showFreq)){
                g <- g + geom_text( aes(label=round(Val,2),y = (max(Val)/2)))
              }
            }  
          } else {
            df  <- df <- BarplotData() %>%
                           group_by(Partition1) %>% 
                           mutate(Freq = Val / sum(Val))
            
            df$Partition2 <- factor(df$Partition2,
                                    levels=rev(levels(df$Partition2)))
            g  <- ggplot(df,
                   aes(x= Partition1,
                       fill=Partition2,
                       y= Freq)) +
              geom_bar( position = "fill",
                        stat = 'identity') + 
              coord_flip() +
              xlab(input$partitionType1) + 
              ylab("Proportion") + 
              scale_fill_manual(values=OrderPartReact()$colPart) + 
              labs(fill = input$partitionType2)
          }
          
        }
      
      g
      # if(input$freq & (input$partitionType2 == "None" | input$wrap)){
      #   g <- g  + geom_text(stat='prop', aes(label=..prop..), position = position_stack(vjust = 0.85))
      #    } else{
      #      g <- g  + geom_text(stat='count', aes(label=..count..), position = position_stack(vjust = 0.85))
      #    }
      #  }
      # 
      # g
      
    })
    
    output$plot_cluster <- renderPlot({
      req(!is.null(Barplot_cluster()))
      Barplot_cluster()
    })
    
    #### Matrix  -----
    #### Metric ----
    
    output$randValue <- renderText({
      req(input$partitionType1)  
      req(input$partitionType2 != "None")
      #req(input$barplot_matrix == "matrix")
      rand <- bluster::pairwiseRand(colData(sce)[,input$partitionType1],
                   colData(sce)[,input$partitionType2],
                   mode="index")
      paste('<b style="color:black;">Rand Index:',round(rand,2),'<b><br>')
    })
    
    PlotMatrix <- reactive({
      req(input$partitionType1)  
      req(input$partitionType2 != "None")
      #req(input$barplot_matrix == "matrix")
      
      if(input$Metric == "Count"){
        tab <- table(colData(sce)[,input$partitionType1],
                     colData(sce)[,input$partitionType2])
        # tab <- log10(tab+10)
        col <- viridis::viridis(100)
      } else if (input$Metric == "Jaccard"){
        tab <- bluster::linkClustersMatrix(colData(sce)[,input$partitionType1],
                                  colData(sce)[,input$partitionType2],
                                  denominator = "union")
        col <- viridis::viridis(100)
      } 
      # else if (input$Metric == "Rand"){
      #   tab <- bluster::pairwiseRand(colData(sce)[,input$partitionType1],
      #                       colData(sce)[,input$partitionType2],
      #                       mode="ratio")
      #   col <- viridis::magma(100) 
      # }
      
      Heatmap(tab, name = input$Metric,
              col = col,
              cluster_rows = F,
              cluster_columns = F,
              show_row_names = T,
              show_column_names = T,
              row_title = input$partitionType1,
              column_title = input$partitionType2,
              cell_fun = if(!(input$showFreq)){set_val(tab = tab)} else{NULL}
      )
      
      
    })
    
    output$plot_matrix <- renderPlot({
    req(!is.null(PlotMatrix()))
      PlotMatrix()
    })
    
    ### Downloads -----
    
    output$export_barplot = downloadHandler(
      filename = function() {"BarPlots_Categories.pdf"},
      content = function(file) {
        pdf(file,
            width = input$pdf_width_barplot,
            height = input$pdf_height_barplot
        )
        Barplot_cluster() %>% plot()
        dev.off()
      }
    )
    
    output$export_matrix = downloadHandler(
      filename = function() {"Matrix_Categories.pdf"},
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
