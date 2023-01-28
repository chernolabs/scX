#####          TOOLS               ######## 
#####        Violin Tab            ######## 
#####                              ######## 

##### ViolinGL Server Module ----
VT_UI <- function(id) {
  tagList(
    fluidRow(
      column(3,
        box(title = htmltools::span(icon("fa-light fa-gears"), " Settings"),
            width = NULL, status = "primary", solidHeader = T, collapsible = F,
          fluidRow(
            column(6,
                   style='padding-left:12px; padding-right:3px;',
                   align="center",
              pickerInput(NS(id,"partitionType"),
                          "Partition",
                          choices = NULL)
            ),
            column(6,
                   style='padding-left:3px; padding-right:12px;',
                   align="center",
              pickerInput(NS(id,"partitionType_wrap"),
                          "Partition wrap",
                          choices = NULL)
            )
          ),
          fileInput(NS(id,'listGenes'),
                    label = NULL,
                    multiple = T, 
                    accept = c("txt/csv", "text/comma-separated-values,
                             text/plain", ".csv", ".xlsx"),
                    buttonLabel = "Search",
                    placeholder = "Select Gene list"),
          htmlOutput(NS(id,"missingGenes")),
          fluidRow(
            column(12,
                   align = "right",
                   style='padding-left:12px; padding-right:12px;',
              uiOutput(NS(id,"previewButton"))
            )
          )
        )
      ),
      column(9,
        box(title = "Preview", width = NULL, solidHeader = TRUE,
            status = "primary",collapsible = T,
          plotOutput(NS(id,"plot")) %>% withSpinner()
        )
      )
    )
  )
}

##### ViolinGL Server Module ----
VT_Server <- function(id,sce) {
  moduleServer(id, function(input, output, session) {
    
    ### Observe Events ----
    updatePickerInput(session,inputId = "partitionType", 
                        choices = names(colData(sce))[sapply(colData(sce), is.factor)])
      
    observeEvent(input$partitionType,ignoreInit = T,{
      req(input$partitionType)
      prt <- names(colData(sce))[sapply(colData(sce), is.factor)]
      updatePickerInput(session,inputId = "partitionType_wrap", 
                        choices = c("None",prt[-(match(input$partitionType,prt))])
      )
    })
    
    ### Gen selected & data preparation ----
    #### Gene List Uploaded  -----
    genes.GL <- reactive({
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
    OrderPartReact <- eventReactive(input$partitionType,{
      req(input$partitionType)
      Col.and.Order(partition = input$partitionType, sce=sce)
    })
    
    ### Preview ------
    output$plot <- renderPlot({
      if(input$partitionType_wrap != "None"){
        df <-data.frame(X=colData(sce)[,input$partitionType],
                        Wrap=colData(sce)[,input$partitionType_wrap])
        g <-ggplot(df, aes(x=X))+
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none") + 
          facet_wrap(.~Wrap) 
      }
      
      else {
        df <-data.frame(X=colData(sce)[,input$partitionType])
        g <- ggplot(df,aes(x=X)) + 
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none") + ggtitle("Preview")
        
      }
      g +
        ylab("log(counts)") +
        xlab(input$partitionType)
      
    })
    
    ### Download ----
    output$previewButton <- renderUI({
      req(length(genes.GL()$genes)>0)
      downloadButton(NS(id,'export'))
    })

    output$export = downloadHandler(
        filename = function() {"Violin_plots.pdf"},
        content = function(file) {
          pdf(file)
          for (i in 1:length(genes.GL()$genes)){
            if(input$partitionType_wrap != "None"){
              df <-data.frame(Y=assay(sce,"logcounts")[genes.GL()$genes[i],],
                              X=colData(sce)[,input$partitionType],
                              Wrap=colData(sce)[,input$partitionType_wrap])
              
              summ <- df  %>% group_by(Wrap,X) %>% summarise(n=n(),Ymax = (max(Y)+0.5))
              
              g <-ggplot(df, aes(y=Y,x=X,fill=X)) + geom_violin(scale="width") + 
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none") + 
                facet_wrap(.~Wrap) + geom_text(aes(label = n,y=Ymax), data = summ) + 
                ggtitle(paste0(genes.GL()$genes[i]," Expression")) +
                scale_fill_manual(values=OrderPartReact()$colPart) +
                ylab("log(counts)") +
                xlab(input$partitionType)
            }
            else {
              df <-data.frame(Y=assay(sce,"logcounts")[genes.GL()$genes[i],],
                              X=colData(sce)[,input$partitionType])
              summ <- df  %>% group_by(X) %>% summarise(n=n(),Ymax = (max(Y)+0.5))
              g <-ggplot(df,aes(y=Y,x=X,fill=X)) + 
                geom_violin(scale="width") + 
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none") + 
                geom_text(aes(label = n, y=Ymax), data = summ) +
                scale_fill_manual(values=OrderPartReact()$colPart) +
                ggtitle(paste0(genes.GL()$genes[i]," Expression")) + 
                ylab("log(counts)") +
                xlab(input$partitionType)
              
            }
            plot(g)
          }
          dev.off()
      }
    )
  
  })
}