#####          TOOLS               ######## 
#####     MultiPlot Tab            ######## 
#####                              ######## 

##### MultiPlot UI Module ----
MultiPlotsUI <- function(id) {
  tagList(
    fluidRow(
      column(3,
        box(title = htmltools::span(icon("gears"), " Settings"),
            width = NULL, status = "primary",solidHeader = T,collapsible = F,
          fluidRow(
            column(6,style='padding-left:12px; padding-right:3px;', align="center",
              pickerInput(NS(id,"plotType"),
                          "  Plot Type",
                          choices = NULL,
                          width = NULL)
            ),
            column(6,style='padding-left:3px; padding-right:12px;', align="center",
              conditionalPanel("input.gene_cluster == 'gene'",ns=NS(id),
                pickerInput(NS(id,"colPal"),
                            "  Col Palette",
                            choices = c("viridis","red","blue"),
                            width = NULL)
              ),
              conditionalPanel("input.gene_cluster == 'cluster'",ns=NS(id),
                pickerInput(NS(id,"colPal_cluster"),
                            "  Col Palette",
                            choices = c("red","blue"),
                            width = NULL)
              )
            )
          ),
          conditionalPanel("input.gene_cluster == 'gene'",ns=NS(id),
            fluidRow(
              column(12,style='padding-left:12px; padding-right:12px;',
                     align="center",
                switchInput(NS(id,"GL_T"),
                            label = "Upload GeneList",
                            size = "small",
                            width = NULL,
                            labelWidth = "100px"),
              )
            ),
            tabsetPanel(id = NS(id,"switcher"), type = "hidden", selected = "panel1",
              tabPanelBody("panel1",
                fluidRow(
                  column(11, style='padding-left:0px; padding-right:2px;',
                   selectizeInput(NS(id,"gen_exp"),
                                   label="Genes",
                                   choices = NULL, 
                                   options = list(maxItems = 10,
                                                  maxOptions = 20,
                                                  placeholder = 'Please select genes to plot'),
                                   width = NULL,
                                   multiple=T)
                  ),
                  column(1, style='padding-left:2px; padding-right:2px; padding-top:4px',
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
                  column(12, style='padding-left:0px; padding-right:0px;',
                    fileInput(NS(id,'listGenes'), label = NULL, multiple = T, 
                              accept = c("txt/csv", "text/comma-separated-values,
                                        text/plain", ".csv", ".xlsx"),
                              buttonLabel = "Search", 
                              placeholder = "Select Gene list"),
                    htmlOutput(NS(id,"missingGenes"))
                  )
                )
              )
            )
          ),
          conditionalPanel("input.gene_cluster == 'cluster'",ns=NS(id),
            fluidRow(
              column(12,style='padding-left:12px; padding-right:12px;',
                pickerInput(inputId = NS(id,"partitionType"),
                            label = "Partition",
                            choices = NULL
                )
              )
            ),
            fluidRow(
              column(10, style='padding-left:12px; padding-right:2px;',
                     selectizeInput(NS(id,"clusterType"), "Clusters",
                                    choices = NULL, 
                                    options = list(placeholder = 'Please select clusters to plot'),
                                    multiple=T)
              ),
              column(2, style='padding-left:2px; padding-right:12px; padding-top:28px',
                     align="center",
                     actionBttn((NS(id,"actionButton")),
                                label = NULL,
                                style = "unite",
                                color = "primary",
                                size = "xs",
                                icon = icon("play"))
              )
            ),
            hr(style = "border-top: 1px solid #0073b7;"),
            fluidRow(align="center",
                     div(style = "display:inline-block;color:#0073b7", icon("circle-info"),
                         title = "To color all clusters simultaneously,\ngo to Cluster Markers tab. The plot toolbar\nincludes the option to download as a png.")
            )
          ),
          conditionalPanel("input.gene_cluster == 'field'",ns=NS(id),
            fluidRow(
              column(10, style='padding-left:12px; padding-right:2px;',
                selectizeInput(NS(id,"fieldType"), "Fields",
                               choices = NULL, 
                               options = list(placeholder = 'Please select fields to plot'),
                               multiple=T)
              ),
              column(2, style='padding-left:2px; padding-right:12px; padding-top:28px',
                     align="center",
                actionBttn((NS(id,"field_actionButton")),
                           label = NULL,
                           style = "unite",
                           color = "primary",
                           size = "xs",
                           icon = icon("play"))
              )
            )
          )
        )
      ),
      column(9,
        tabBox(id = NS(id,"gene_cluster"),
               selected = "gene",
               width = NULL,
          tabPanel("by Gene Expression",value = "gene",
                   dropdownButton(
					fluidRow(
						column(4, style='padding-left:12px; padding-right:3px;', numericInput(NS(id,"pdf_width_bygene"),"Width",value = 7)),
						column(4, style='padding-left:3px; padding-right:3px;', numericInput(NS(id,"pdf_height_bygene"),"Height",value = 7)),
						column(4, style='padding-left:3px; padding-right:12px; padding:16px', downloadButton(NS(id,'export_bygene')))
					),
                     circle = FALSE,
                     status = "primary",
                     icon = icon("fa-thin fa-download"),
                     width = "300px",
                     size= "sm",
                     up = F,
                     tooltip = tooltipOptions(title = "Press to Download")
                   ),
            plotOutput(NS(id,"plot_gene"),
                       height = "100vh") %>% withLoader(type='html',loader = 'dnaspin')
          ),
          tabPanel("by Cluster",value = "cluster",
                   dropdownButton(
					fluidRow(
						column(4, style='padding-left:12px; padding-right:3px;', numericInput(NS(id,"pdf_width_bycluster"),"Width",value = 7)),
						column(4, style='padding-left:3px; padding-right:3px;', numericInput(NS(id,"pdf_height_bycluster"),"Height",value = 7)),
						column(4, style='padding-left:3px; padding-right:12px; padding:16px', downloadButton(NS(id,'export_bycluster')))
					),
                     circle = FALSE,
                     status = "primary",
                     icon = icon("fa-thin fa-download"),
                     width = "300px",
                     size= "sm",
                     up = F,
                     tooltip = tooltipOptions(title = "Press to Download")
                   ),
            plotOutput(NS(id,"plot_cluster"),
                      height = "100vh") %>% withLoader(type='html',loader = 'dnaspin')
          ),
          tabPanel("by Field",value = "field",
                   dropdownButton(
					fluidRow(
						column(4, style='padding-left:12px; padding-right:3px;', numericInput(NS(id,"pdf_width_byfield"),"Width",value = 7)),
						column(4, style='padding-left:3px; padding-right:3px;', numericInput(NS(id,"pdf_height_byfield"),"Height",value = 7)),
						column(4, style='padding-left:3px; padding-right:12px; padding:16px', downloadButton(NS(id,'export_byfield')))
					),
                     circle = FALSE,
                     status = "primary",
                     icon = icon("fa-thin fa-download"),
                     width = "300px",
                     size= "sm",
                     up = F,
                     tooltip = tooltipOptions(title = "Press to Download")
                   ),
                   plotOutput(NS(id,"plot_field"),
                              height = "100vh") %>% withLoader(type='html',loader = 'dnaspin')
          )
        )
      )
    )
  )
}

##### MultiPlot Server Module ----
MultiPlotsServer <- function(id,sce) {
  moduleServer(id, function(input,output,session) {
    
    ### Observe Events ----
    observeEvent(ignoreInit = T,input$GL_T,{
      if(input$GL_T) {
        updateTabsetPanel(inputId = "switcher", selected = "panel2")
      } else{
        updateTabsetPanel(inputId = "switcher", selected = "panel1")
      }
    })
    
    updateSelectizeInput(session, 'gen_exp', choices = rownames(sce), server = TRUE)
 
    dimVector <- reactive({
      sapply(reducedDims(x = sce),FUN = ncol) 
    })
    
    observeEvent(dimVector(),{
      updatePickerInput(session,inputId = "plotType", choices = rev(names(which(dimVector() == 2  | dimVector() > 3))))
    })
    
    updatePickerInput(session, 'partitionType', 
                      choices = names(colData(sce))[sapply(colData(sce), is.factor)])
    
    observeEvent(input$partitionType,{
      updateSelectizeInput(session,inputId = "clusterType",
                           choices = levels(colData(sce)[,input$partitionType]),
                           server = TRUE)
    }) 
    
    updateSelectizeInput(session, 'fieldType', 
                         choices = names(colData(sce))[sapply(colData(sce), is.numeric)],
                         server = TRUE)
    
    ### Gen selected & data preparation ----
    #### Selected Gene -----
    
    genes.L <- eventReactive(input$action,{
      req(input$GL_T == F)
      req(input$gen_exp)
      input$gen_exp
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
    
    
    #### Plots ----
    list_plots <- reactive({
      req(input$plotType)
      if(input$GL_T){
        req(length(genes.GL()$genes) > 0)
        feature <- genes.GL()$genes
      } else { 
        req(!is.null(genes.L()))
        feature <- genes.L() 
      }
      
       list_plots <- list()
       for(i in 1:length(feature)){
         list_plots[[feature[i]]] <- scater::plotReducedDim(object = sce, dimred = input$plotType,ncomponents = 2, colour_by=feature[i]) + ggtitle(feature[i]) 
         
       }
      list_plots
    })
    
    list_gene_plot <- reactive({
      f_plot <- list_plots()
      if(input$colPal != "viridis"){
      for(i in 1:length(f_plot)){
        f_plot[[i]] <- f_plot[[i]] +  scale_color_continuous(low="lightgrey",high=input$colPal) + labs(color=names(f_plot)[i])
      }
      }
      
      cowplot::plot_grid(plotlist = f_plot)
    })
      
    output$plot_gene <- renderPlot({
      req(!is.null(list_gene_plot()))
      list_gene_plot()
      })
    
    ### Clusters ----
    
    list_cluster_plot <- eventReactive(c(input$actionButton,input$plotType,input$colPal_cluster),{
      req(length(input$clusterType) > 0)
      req(!is.null(input$plotType))
      f_plot <- list()
      for(i in 1:length(input$clusterType)){
        f_plot[[i]] <- reducedDimPlot_cluster(sce = sce,reducedDim = input$plotType,partition = input$partitionType,
                                              cluster = input$clusterType[i],alpha = 0.5,palette = input$colPal_cluster)
      }
      cowplot::plot_grid(plotlist = f_plot)
      
    })
    
    output$plot_cluster <- renderPlot({
      req(!is.null(list_cluster_plot()))
      list_cluster_plot()
    })
    
    ### Fields ---
    
    list_field_plot <- eventReactive(c(input$field_actionButton,input$plotType),{
      req(length(input$fieldType) > 0)
      req(!is.null(input$plotType))
      list_plots <- list()
      feature <- input$fieldType
      for(i in 1:length(feature)){
        list_plots[[feature[i]]] <- scater::plotReducedDim(object = sce, dimred = input$plotType,ncomponents = 2, colour_by=feature[i]) + ggtitle(feature[i]) + 
          scale_color_distiller(palette ='YlOrRd',direction = 1) + labs(color=feature[i])
        
      }
      
      cowplot::plot_grid(plotlist = list_plots)
    })
    
    
    output$plot_field <- renderPlot({
      req(!is.null(list_field_plot()))
      list_field_plot()
    })
    
  
    ### Downloads -----
    
    output$export_bygene = downloadHandler(
      filename = function() {"MultiPlot_byGeneExpression.pdf"},
      content = function(file) {
        pdf(file,
            width = input$pdf_width_bygene,
            height = input$pdf_height_bygene
        )
        list_gene_plot() %>% plot()
        dev.off()
      })
    
    output$export_bycluster = downloadHandler(
      filename = function() {"MultiPlot_byCluster.pdf"},
      content = function(file) {
        pdf(file,
            width = input$pdf_width_bycluster,
            height = input$pdf_height_bycluster
        )
        list_cluster_plot() %>% plot()
        dev.off()
      })
    
    output$export_byfield = downloadHandler(
      filename = function() {"MultiPlot_byField.pdf"},
      content = function(file) {
        pdf(file,
            width = input$pdf_width_byfield,
            height = input$pdf_height_byfield
        )
        list_field_plot() %>% plot()
        dev.off()
      })
    
  })
}

