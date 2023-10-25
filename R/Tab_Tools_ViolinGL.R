#####          TOOLS               ######## 
#####        Violin Tab            ######## 
#####                              ######## 

##### ViolinGL Server Module ----
VT_UI <- function(id) {
  tagList(
    fluidRow(
      column(3,
        box(title = htmltools::span(icon("gears"), " Settings"),
            width = NULL, status = "primary", solidHeader = T, collapsible = F,
          fluidRow(
            column(6,
                   style='padding-left:12px; padding-right:3px;',
                   align="center",
              pickerInput(NS(id,"partitionType"),
                          "Category",
                          choices = NULL)
            ),
            column(6,
                   style='padding-left:3px; padding-right:12px;',
                   align="center",
              pickerInput(NS(id,"partitionType_wrap"),
                          "Category wrap",
                          choices = NULL)
            )
          ),
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
                                 label=NULL,
                                 choices = NULL, 
                                 options = list(maxItems = 10,
                                 maxOptions = 20,
                                 placeholder = 'Please select genes to plot'),
                                 width = NULL,
                                 multiple=T
                  )
                ),
                column(1, style='padding-left:2px; padding-right:2px; padding-top:4px',
                  actionBttn(NS(id,"action"),
                             label = NULL,
                             style = "unite",
                             color = "primary",
                             size = "xs",
                             icon = icon("play")
                  )
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
          ),
          uiOutput(NS(id,"previewButton"))
        )
      ),
      column(9,
        box(title = "Preview - Downloaded plots will have the following layout",
            width = NULL,
            solidHeader = T,
            collapsible = F,
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
    
    observeEvent(ignoreInit = T,input$GL_T,{
      if(input$GL_T) {
        updateTabsetPanel(inputId = "switcher", selected = "panel2")
      } else{
        updateTabsetPanel(inputId = "switcher", selected = "panel1")
      }
    })
    
    updateSelectizeInput(session, 'gen_exp', choices = rownames(sce), server = TRUE)
    
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
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
        
      }
      g +
        ylab("log(counts)") +
        xlab(input$partitionType)
      
    })
    
    ### Download ----
    output$previewButton <- renderUI({
      if(input$GL_T){
        req(length(genes.GL()$genes) > 0)
      } else { 
        req(!is.null(genes.L()))
      }
      tagList(
        hr(),
        fluidRow(
          column(3,style='padding-left:9px; padding-right:1px;',align="center",
                 numericInput(NS(id,"pdf_width"),"Width",value = 7)),
          column(3,style='padding-left:1px; padding-right:3px;',align="center",
                 numericInput(NS(id,"pdf_height"),"Height",value = 7)),
          column(6,style='padding-left:3px; padding-right:12px;padding-top:25px;',align="center",
                 downloadButton(NS(id,'export'))
          ) 
        )
      )
    })

    output$export = downloadHandler(
        filename = function() {"Violin_plots.pdf"},
        content = function(file) {
          if(input$GL_T){
            req(length(genes.GL()$genes) > 0)
            feature <- genes.GL()$genes
          } else { 
            req(!is.null(genes.L()))
            feature <- genes.L() 
          }
          pdf(file,
              width = input$pdf_width,
              height = input$pdf_height
              )
          for (i in 1:length(feature)){
            if(input$partitionType_wrap != "None"){
              df <-data.frame(Y=assay(sce,"logcounts")[feature[i],],
                              X=colData(sce)[,input$partitionType],
                              Wrap=colData(sce)[,input$partitionType_wrap])
              
              summ <- df  %>% group_by(Wrap,X) %>% summarise(n=n(),Ymax = (max(Y)+0.5))
              
              g <-ggplot(df, aes(y=Y,x=X,fill=X)) + geom_violin(scale="width") + 
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none") + 
                facet_wrap(.~Wrap) + geom_text(aes(label = n,y=Ymax), data = summ) + 
                ggtitle(paste0(feature[i]," Expression")) +
                scale_fill_manual(values=OrderPartReact()$colPart) +
                ylab("log(counts)") +
                xlab(input$partitionType)
            }
            else {
              df <-data.frame(Y=assay(sce,"logcounts")[feature[i],],
                              X=colData(sce)[,input$partitionType])
              summ <- df  %>% group_by(X) %>% summarise(n=n(),Ymax = (max(Y)+0.5))
              g <-ggplot(df,aes(y=Y,x=X,fill=X)) + 
                geom_violin(scale="width") + 
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none") + 
                geom_text(aes(label = n, y=Ymax), data = summ) +
                scale_fill_manual(values=OrderPartReact()$colPart) +
                ggtitle(paste0(feature[i]," Expression")) + 
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
