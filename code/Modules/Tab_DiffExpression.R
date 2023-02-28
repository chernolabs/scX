#####                              ######## 
#####  DIFFERENTIAL EXPRESSION Tab ######## 
#####                              ######## 

### Volcano UI Module -----
VolcanoUI <- function(id) {
  tagList(
    fluidRow(
      #### Side panel -----
      column(width = 4,
             box(
               title = htmltools::span(icon("fa-light fa-gears"), " Settings"), width = NULL,
               status = "primary",solidHeader = T,collapsible = T,
               pickerInput(
                 inputId = NS(id,"partitionType"),
                 label = "Partition", 
                 choices = NULL
               ),
               selectizeInput(NS(id,"clusterType"), "Clusters",
                              choices = NULL, 
                              options = list(maxItems = 2,
                                             placeholder = 'Please select two options below'),
                              multiple=T),
               hr(style = "border-top: 1px solid #0073b7;"),
               fluidRow(
                 column(4,style='padding-left:12px; padding-right:5px;',align="center",
                        numericInput(NS(id,"fdr_thr"),"FDR",value = 0.05)),
                 column(4,style='padding-left:5px; padding-right:5px;',align="center",
                        numericInput(NS(id,"fc_thr"),"Fold Change",value = 2)),
                 column(4,style='padding-left:5px; padding-right:12px;padding-top:25px;',align="center",
                        actionBttn((NS(id,"button")),label = "Calculate", style = "unite",color = "primary",size = "sm"))
               )
             ),
             conditionalPanel("typeof output.VolcanoPlot !== 'undefined'", ns = NS(id),
                              box(title = NULL, width = NULL,solidHeader = F,
                                  checkboxGroupButtons(inputId = NS(id,"box_exp"), 
                                                       label = NULL,
                                                       choices = c("DOWN","UP","Selected"),
                                                       selected = c("DOWN","UP"),
                                                       direction = "horizontal",
                                                       justified = TRUE,
                                                       checkIcon = list(yes = icon("ok", 
                                                                                   lib = "glyphicon"))
                                  ),
                                  DTOutput(NS(id,"DTVolcano")) %>% withLoader(type='html',loader = 'loader6')
                              )
             )
      ),
      ## Main Panel -----
      column(width = 8,
             box(title = "Volcano Plot", width = NULL, solidHeader = TRUE, collapsible = T,
                 plotlyOutput(NS(id,"VolcanoPlot"),height = "100vh") %>% withLoader(type='html',loader = 'dnaspin')
             ),
             conditionalPanel("typeof output.VolcanoPlot !== 'undefined'", ns = NS(id),
                              box(title="Expression Plots",width = NULL, solidHeader = T, collapsible=T,
                                  uiOutput(NS(id,"Violin.Bar_Input")),
                                  conditionalPanel("input.Cell_Exp == 'Heatmap'", ns = NS(id),
                                                   dropdownButton(
                                                     prettySwitch(NS(id,"norm_heat"),"Norm per gene",value = F,status = "success",
                                                                  fill = TRUE),
                                                     prettySwitch(NS(id,"mean_heat"),"Mean per group",value = F,status = "success",
                                                                  fill = TRUE),
                                                     prettySwitch(NS(id,"cluster_row"),"Cluster Row",value = F,status = "success",
                                                                  fill = TRUE),
                                                     prettySwitch(NS(id,"cluster_column"),"Cluster Column",value = F,status = "success",
                                                                  fill = TRUE),
                                                     circle = FALSE, status = "primary", icon = icon("cog"), width = "300px", size= "sm",
                                                     tooltip = tooltipOptions(title = "Press to see Heatmaps setting")
                                                   )
                                  ),
                                  conditionalPanel("input.Cell_Exp == 'Dotplot'", ns = NS(id),
                                                   dropdownButton(
                                                     prettySwitch(NS(id,"ord_dotplot"),"Cluster Row",value = F,status = "success",
                                                                  fill = TRUE),
                                                     prettySwitch(NS(id,"scale_dotplot"),"Scale",value = F,status = "success",
                                                                  fill = TRUE),
                                                     prettySwitch(NS(id,"center_dotplot"),"Center",value = F,status = "success",
                                                                  fill = TRUE),
                                                     circle = FALSE, status = "primary", icon = icon("cog"), width = "300px",
                                                     size= "sm",
                                                     tooltip = tooltipOptions(title = "Press to see DotPlot setting")
                                                   )
                                  ), 
                                  uiOutput(NS(id,"Cell_Plots")) %>% withSpinner()
                              ) 
             )  
      )
    )
  )
}



### Volcano Server Module -----
VolcanoServer <- function(id,sce,sce.markers) {
  moduleServer(id, function(input,output,session) {
    
    ### Observe Events ----
    updatePickerInput(session, 'partitionType', 
                      choices = names(sce.markers))
    
    observeEvent(input$partitionType,{
      updateSelectizeInput(session,inputId = "clusterType",choices = names(sce.markers[[input$partitionType]]),server = TRUE)
    }) 
    
    ### Volcano Plot ----
    #Preparing data 
    VolcanoC <- eventReactive(c(input$button,input$partitionType),{
      req(input$partitionType)
      req(length(input$clusterType) > 1)
      
      dta <- sce.markers[[input$partitionType]][[input$clusterType[2]]][[c(paste0("stats.",input$clusterType[1]))]]
      dta$diffexpressed <- "NO"
      fc_Thr <- abs(log2(input$fc_thr))
      fdr_Thr <- abs(log(input$fdr_thr))
      dta$diffexpressed[dta$logFC > fc_Thr & -dta$log.FDR > fdr_Thr] <- "UP"
      dta$diffexpressed[dta$logFC < -fc_Thr & -dta$log.FDR > fdr_Thr] <- "DOWN"
      
      list(Exp=as.data.frame(dta),
           Clusters = input$clusterType,
           Thr_fc = fc_Thr,
           Thr_fdr = fdr_Thr,
           UP = rownames(dta)[dta$diffexpressed == "UP"],
           DOWN = rownames(dta)[dta$diffexpressed == "DOWN"]
      )
    })
    
    #Plotting
    output$VolcanoPlot <- renderPlotly({
      req(!is.null(VolcanoC()))
      
      plot_ly(source = "VolcanoPlot",
              data = VolcanoC()$Exp, x = ~logFC, y = ~-log.FDR, 
              color= ~ diffexpressed, colors=c("UP"="red","NO"= "black","DOWN"= "blue"),
              customdata = ~I(row.names(VolcanoC()$Exp)),
              hovertemplate = paste('<b>%{text}</b>',
                                    '<br> log2(FC): %{x:.2f}',
                                    '<br>-log(FDR): %{y:.2f}'
              ),
              text = ~rownames(VolcanoC()$Exp)) %>%  
        layout(title = paste(VolcanoC()$Clusters[1],"→",VolcanoC()$Clusters[2]),
               xaxis = list(title = 'log2(FC)'),
               yaxis = list(title = '-log(FDR)'),
               shapes = list(vline(x = VolcanoC()$Thr_fc),
                             vline(x = -VolcanoC()$Thr_fc),
                             hline(y = VolcanoC()$Thr_fdr)
               )) %>% 
        event_register("plotly_selecting") %>% 
        toWebGL()
    })
    
    #DT ----
    output$DTVolcano <- renderDT(server = FALSE,
                                 datatable({
                                   VolcanoC()$Exp[Genes(),c("diffexpressed","logFC","log.FDR")]
                                 },
                                 extensions = 'Buttons',
                                 options = list(language = list(zeroRecords = "No genes selected"),
                                                dom = 'Bfrtip',
                                                exportOptions = list(header = ""),
                                                buttons = c('copy', 'csv', 'excel', 'pdf'),
                                                rowCallback = DT::JS(js_volcano),
                                                scrollX=T,
                                                autoWidth = F),
                                 rownames= TRUE,
                                 selection = 'none',
                                 caption = htmltools::tags$caption(
                                   style = paste0('caption-side: top; text-align: center; font-weight: bold;color:black;background-image: linear-gradient(to right,',
                                                  OrderPartReact()$colPart[[VolcanoC()$Clusters[1]]],',',
                                                  OrderPartReact()$colPart[[VolcanoC()$Clusters[2]]],');'
                                   ),
                                   paste0(VolcanoC()$Clusters[1], "→",VolcanoC()$Clusters[2]))
                                 ) %>% 
                                   formatStyle('diffexpressed',
                                               color = styleEqual(c("UP", "NO","DOWN"), c("red", "black", "blue"))
                                   )
    )
    
    #Selecting genes from DT ----
    Selected <- reactive({
      d <- event_data(source = "VolcanoPlot","plotly_selected")
      d <- d[,"customdata"]
      if (is.null(d)) NULL else d
    })
    
    Genes <- reactive({
      req(!is.null(VolcanoC()))
      req(input$box_exp)
      UP <- DOWN <- Select <- c()
      if("UP" %in% input$box_exp) {UP <- VolcanoC()$UP} 
      if("DOWN" %in% input$box_exp) {DOWN <- VolcanoC()$DOWN} 
      if("Selected" %in% input$box_exp) {Select <- Selected()} 
      
      genes <- unique(c(UP,DOWN,Select))
      genes
    })
    
    
    ##### Expression Plots ----
    
    #### Volcano&SpikePlots ----
    
    ExpressionL <- eventReactive(Genes(),{
      req(length(Genes()) > 0)
      diff_Exp <- factor(VolcanoC()$Exp[Genes(),"diffexpressed",drop=T], levels =c("UP","DOWN","NO"))
      diff_Exp <- droplevels(diff_Exp)
      if(length(Genes()) > 1) {
        exp_vtor <- apply(logcounts(sce)[Genes(),,drop=F],2,function(x){
          tapply(x,diff_Exp,mean)
        })
        exp_vtor <- t(exp_vtor) %>% melt(measure.vars=rownames(exp_vtor),
                                         value.name="Y")
        colnames(exp_vtor)[1:2] <- c("CellID","DiffExpression")
      } else {
        exp_vtor <- logcounts(sce)[Genes(),,drop=F] %>% as.matrix
        exp_vtor<- t(exp_vtor) %>% as.data.frame
        colnames(exp_vtor) <- "Y"
        exp_vtor$DiffExpression <- diff_Exp
      }
      
      list(Exp = exp_vtor,Genes=Genes(), diff_Exp = diff_Exp)
    })
    
    ### Radio Input ----
    # Colors and Order 
    OrderPartReact <- eventReactive(input$partitionType,{
      req(input$partitionType)
      Col.and.Order(partition = input$partitionType, sce=sce)
    })
    
    output$Violin.Bar_Input <- renderUI({
      req(!is.null(ExpressionL()))
      radioGroupButtons(inputId = NS(id,"Cell_Exp"), label=NULL,
                        choices = c("Violin", "SpikePlot","Heatmap","Dotplot"),
                        direction = "horizontal",justified = T,individual=T)
    })
    
    ViolinReact <- reactive({
      req(!is.null(ExpressionL()) & input$Cell_Exp == "Violin")
      df <- ExpressionL()$Exp
      df$X <- factor(colData(sce)[,input$partitionType])
      df
    })
    
    ViolinReact_Cell <- reactive({
      ViolinReact() %>% group_by(X,DiffExpression) %>% summarise(n=n(),Ymax = (max(Y)+0.5))
    })
    
    ViolinPlot <-reactive({
      req(!is.null(ExpressionL()) & input$Cell_Exp == "Violin")
      ggplot(ViolinReact()) + 
        geom_violin(aes(y = Y, 
                        x = X, 
                        fill = X),
                    scale="width") + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              legend.position = "none") + 
        xlab("Clusters") + 
        ylab("log(counts)") +
        facet_wrap(.~DiffExpression,ncol = 1) +
        ggtitle("Mean Expression") +
        scale_fill_manual(values=OrderPartReact()$colPart) +
        geom_text(aes(label = n,x=X, y=Ymax),
                  data = ViolinReact_Cell()) 
    })
    
    SpikePlot <-reactive({
      req(!is.null(ExpressionL()) & input$Cell_Exp == "SpikePlot")
      df <- ExpressionL()$Exp
      lt <- df  %>% split(df$DiffExpression)
    })
    
    output$Violin.Bar_Plot <- renderPlot({
      req(!is.null(ExpressionL()))
      req(input$Cell_Exp)
      if(input$Cell_Exp == "Violin"){
        ViolinPlot()
      } else { 
        par(mfrow=c((length(SpikePlot())),1))
        for(i in 1:length(SpikePlot())){
          barplot(SpikePlot()[[i]][OrderPartReact()$ordPart,"Y",drop=T],
                                    col = OrderPartReact()$colPart[colData(sce)[,input$partitionType]][OrderPartReact()$ordPart],
                                    border = OrderPartReact()$colPart[colData(sce)[,input$partitionType]][OrderPartReact()$ordPart],
                                    ylab = "log(counts)", main = names(SpikePlot())[i], names.arg = F) 
          abline(h=mean(SpikePlot()[[i]][SpikePlot()[[i]]$Y>0,"Y",drop=T]),lty=2,col="grey")
        }
        legend("bottom", legend = names(OrderPartReact()$colPart), col = OrderPartReact()$colPart,
               pch=19, ncol=6, xpd=T, inset=c(0,-0.10))
        par(mfrow=c(1,1))
      }
      
    })
    
    ####  Heatmap ----
    
    HeatmapL <- eventReactive(c(Genes(),input$norm_heat),{
      req(length(Genes()) > 0)
      req(input$Cell_Exp == "Heatmap")
      value <- ifelse(input$norm_heat,yes = "logcounts.norm",no = "logcounts")
      exp_mtx <-as.matrix(assay(sce,value)[Genes(),])
      if(length(Genes()) == 1) {
        exp_mtx <- t(exp_mtx)
      }
      rownames(exp_mtx) <- Genes()
      diff_Exp <- factor(VolcanoC()$Exp[Genes(),"diffexpressed",drop=T], levels =c("UP","DOWN","NO"))
      color_genes <- setNames(nm = c("UP", "NO","DOWN"),object =  c("red", "black", "blue"))[as.vector(diff_Exp)]
      
      c(list(Exp=exp_mtx,
             color_genes=color_genes,
             diff_exp = diff_Exp),Genes())
    })
    
    Heatmap_Plot <- eventReactive(c(HeatmapL(),input$partitionType,input$mean_heat,input$cluster_row,input$cluster_column),{
      req(input$partitionType)
      req(input$Cell_Exp == "Heatmap")
      req(!is.null(HeatmapL()))
      req(!is.null(ExpressionL()))
      req(!is.null(OrderPartReact()))
      
      if(input$mean_heat) {
        dta <- apply(HeatmapL()$Exp,1,FUN =  function(x){
          tapply(X = x,INDEX = colData(sce)[,input$partitionType],FUN=mean)
        }) %>% t
        ht <-   HeatmapAnnotation(Type = levels(colData(sce)[,input$partitionType]),
                                  col=list(Type=OrderPartReact()$colPart),
                                  annotation_legend_param = list(Type = list(title = input$partitionType)
                                  ),
                                  annotation_label = c(input$partitionType),
                                  show_legend = c(Type =FALSE),
                                  show_annotation_name = T)
        
        h1 <- Heatmap(dta,
                      col = viridis(100),
                      row_names_gp = gpar(col=HeatmapL()$color_genes),
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
                      # row_title = if(nrow(HeatmapL()$Exp) < 50){ "Genes" }else{ F},
                      row_gap = unit(1, "mm"),
                      column_gap = unit(1, "mm"),
                      show_row_names = if(nrow(HeatmapL()$Exp) < 50){ T }else{ F},
                      # show_column_names = T,
                      top_annotation = ht,
                      row_split = HeatmapL()$diff_exp,
                      cluster_row_slices = F,
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
        h1 <- Heatmap(as.matrix(HeatmapL()$Exp),
                      row_names_gp = gpar(col=HeatmapL()$color_genes),
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
                      # row_title = if(nrow(HeatmapL()$Exp) < 50){ "Genes" }else{ F},
                      row_gap = unit(1, "mm"),
                      column_gap = unit(1, "mm"),
                      show_row_names = if(nrow(HeatmapL()$Exp) < 50){ T }else{ F},
                      show_column_names = F,
                      top_annotation = ht,
                      column_split = colData(sce)[,input$partitionType],
                      row_split = HeatmapL()$diff_exp,
                      cluster_column_slices = F,
                      cluster_row_slices = F,
                      use_raster = TRUE,
                      raster_by_magick = TRUE
        )
      }
      h1
    })
    
    output$plot_heatmap <- renderPlot({
      req(input$Cell_Exp == "Heatmap")
      Heatmap_Plot()
    })
    ####  Dotplots ----
    output$plot_DotPlot <- renderPlotly({
      req(input$Cell_Exp == "Dotplot")
      req(input$partitionType)
      req(length(Genes()) > 0)
      feature <- Genes()
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
      diff_Exp <- factor(VolcanoC()$Exp[as.character(g$data$Feature),"diffexpressed",drop=T], levels =c("UP","DOWN","NO"))
      diff_Exp <- droplevels(diff_Exp)
      g$data$Diff <- diff_Exp
      g <- g + facet_wrap(.~Diff,scales="free_y",ncol=1)
      ggplotly(g)
    })
    
    ####  Final output ---- 
    output$Cell_Plots <- renderUI({
      req(input$Cell_Exp)
      if(input$Cell_Exp == "Heatmap"){
        plotOutput(NS(id,"plot_heatmap"),height = "100vh")
      }
      else if(input$Cell_Exp == "Dotplot"){
        plotlyOutput(NS(id,"plot_DotPlot"),height = "100vh")
      }
      else{
        plotOutput(NS(id,"Violin.Bar_Plot"),height = "100vh") 
      }
    })
    
    
  })
}

