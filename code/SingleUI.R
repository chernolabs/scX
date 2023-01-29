### Packages ----
require(shinydashboard)
require(shiny)
require(dplyr)
require(tidyr)
require(ggplot2)
require(highcharter)
require(htmltools)
require(lubridate)
require(stringr)
require(withr)
require(treemap)
require(cowplot)
require(DT)
require(shinyBS)
require(shinyjs)
require(WDI)
require(geosphere)
require(magrittr)
require(shinyWidgets)
require(timevis)
require(bslib)
require(showtext)
require(thematic)
require(shinycssloaders)
require(shinycustomloader)
require(shinydisconnect)
options(spinner.color="#006272")
theme_set(theme_bw())


# nombre_dataset <- "scXplorer"
path_git <- "https://www.leloir.org.ar/biologia-de-sistemas-integrativa?area=bioinformatica-y-biologia-computacional"
path_button <- "https://www.youtube.com/watch?v=OeVv7zc358E"

#### ShinyDashboard ----
  #### Header ----
  header <- 
    dashboardHeader(title = HTML(nombre_dataset), 
                    disable = FALSE, 
                    titleWidth  = 250)
  
  #### SideBar ----
  siderbar <- 
    dashboardSidebar(
                     width = 250,
                     useShinyjs(),
      sidebarMenu(
                  id = 'sidebar',
                  style = "position: relative; overflow: visible;",
                  HTML(
                    paste0(
                      "<br>","<a href='",
                      path_git,
                      "' target='_blank'><img style =
                      'display: block; margin-left: auto; margin-right: auto;'
                      src='scXplorer-03.png' width = '186'></a>",
                      "<br>",
                      paste0(
                        "<p style = 'text-align: center;'><small><a href='",
                        path_button,
                        "' target='_blank'>Click here to listen to a temaiken</a></small></p>"
                      ),
                      "<br>"
                    )
                  ),
                  tags$style("@import url(https://use.fontawesome.com/releases/v6.2.0/css/all.css);"),
                  menuItem("Global Options", tabName = "globalTab", icon = icon("fa-regular fa-globe")),
                    div( id = 'sidebar_global',
                      conditionalPanel("input.sidebar === 'globalTab'",
                        sliderInput(inputId = "point.size",
                                    label = "Point size ScatterPlots",
                                    min = 5,max = 50,value = 20),
                        fluidRow(column=12, align = "right",
                                 style='padding-left:12px; padding-right:12px;',
                          actionBttn(inputId = "button.config",
                                     label = "Apply", 
                                     style = "stretch",
                                     color = "primary")
                        )
                      )
                    ),
                  menuItem("Summary", tabName = "smryTab", icon = icon("fa-regular fa-bookmark"),selected = T),              
      menuItem("Markers", tabName = "markersTab",
               icon = icon('fa-solid fa-location-dot'), startExpanded = F,
        menuSubItem('Cluster markers', tabName = "markers_cluster",
                    icon = icon('fa-solid fa-map-location-dot')
        ),
        menuSubItem('Find new markers', tabName = "markers_new",
                    icon = icon('fa-solid fa-magnifying-glass-location')
        )
      ),
      menuItem("Gene Expression", tabName = "g_expTab",
               icon = icon('fa-solid fa-chart-simple'), startExpanded = F,
        menuSubItem('Expression', tabName = "g_exp",
                    icon = icon('fa-solid fa-signal')
        ),
        menuSubItem('Co-expression', tabName = "g_coexp",
                    icon = icon('fa-regular fa-clone'))
        ),
      menuItem("Differential Expression", tabName = "volcanoTab",
               icon = icon("fa-solid fa-chart-column")
      ),
      menuItem("Partitions", tabName = "clustersTab",
               icon = icon("fa-solid fa-circle-nodes")
      ),
      menuItem("Tools", tabName = "toolsTab",
               icon = icon("fa-solid fa-toolbox"), startExpanded = F,
        menuSubItem('Violin Gene List', tabName = "t_VGL",
                    icon = icon('fa-solid fa-wrench')
        ),
        menuSubItem('MultiPlots', tabName = "t_MP",
                    icon = icon('fa-solid fa-hammer')
        )
      )
    )
  )
  #### Body ----
  body <- 
    dashboardBody( 
      useSweetAlert(),
      
      ### Styling ----
      tags$head(
        tags$link(rel="shortcut icon", href="scXplorer-03.ico"),
        tags$script(glue('document.title = {nombre_dataset}')),
        tags$style(
          HTML(".tab-content { padding-left: 20px; padding-right: 30px; }")
        ),
        tags$style(
          HTML('/* change size of icons in sub-menu items */
               .sidebar .sidebar-menu .treeview-menu>li>a>.fa {
               font-size: 15px;
               }
               .sidebar .sidebar-menu .treeview-menu>li>a>.glyphicon {
               font-size: 13px;
               }
               /* Hide icons in sub-menu items */
               .sidebar .sidebar-menu .treeview>a>.fa-angle-left {
               display: none;
               }'
          )
        ),
        tags$style(
          HTML("hr {border-top: 1px solid #000000;}")
        ),
        ## to not show error message in shiny
        tags$style(
          HTML(".shiny-output-error { visibility: hidden; }")
        ),
        tags$style(
          HTML(".shiny-output-error:before { visibility: hidden; }")
        ),
        ## heand dropdown menu size
        tags$style(
          HTML('.navbar-custom-menu>.navbar-nav>li:last-child>.dropdown-menu 
               { width:10px; font-size:10px; padding:1px; margin:1px;}')
        ),
        tags$style(HTML('.navbar-custom-menu> .navbar-nav> li:last-child > .dropdown-menu > h4 
                        {width:0px; font-size:0px; padding:0px; margin:0px;}')
        ),
        tags$style(HTML('.navbar-custom-menu> .navbar-nav> li:last-child > .dropdown-menu > p 
                        {width:0px; font-size:0px; padding:0px; margin:0px;}')
        )
      ),
      ### Main -----
      tabItems(
        tabItem(tabName = "smryTab",
          QC_UI(id="qc")
        ),
        tabItem(tabName = "markers_cluster",
          HTML('<h4><span style="font-family:Trebuchet MS,Helvetica,sans-serif"><strong><span style="color:#4e5f70">Cluster markers</span></strong></span></h4>'),
          markersUI(id="markers")
        ),
        tabItem(tabName = "markers_new",
          HTML('<h4><span style="font-family:Trebuchet MS,Helvetica,sans-serif"><strong><span style="color:#4e5f70">Find new markers</span></strong></span></h4>'),
          N_markersUI(id="n_markers")
        ),
        tabItem(tabName = "g_exp",
          HTML('<h4><span style="font-family:Trebuchet MS,Helvetica,sans-serif"><strong><span style="color:#4e5f70">Expression</span></strong></span></h4>'),
          ExpressionUI(id="Exp")
        ),
        tabItem(tabName = "g_coexp",
          HTML('<h4><span style="font-family:Trebuchet MS,Helvetica,sans-serif"><strong><span style="color:#4e5f70">Co-Expression</span></strong></span></h4>'),
          COExpUI(id="Co-exp")
        ),
        tabItem(tabName = "volcanoTab",
          HTML('<h4><span style="font-family:Trebuchet MS,Helvetica,sans-serif"><strong><span style="color:#4e5f70">Differential Expression</span></strong></span></h4>'),
          VolcanoUI(id="volcano")
        ),
        tabItem(tabName = "clustersTab",
          HTML('<h4><span style="font-family:Trebuchet MS,Helvetica,sans-serif"><strong><span style="color:#4e5f70">Partition Analysis</span></strong></span></h4>'),
          Clusters_UI("cluster")
        ),
        tabItem(tabName = "t_VGL",
           HTML('<h4><span style="font-family:Trebuchet MS,Helvetica,sans-serif"><strong><span style="color:#4e5f70">Violin by Gene list </span></strong></span></h4>'),
           VT_UI(id = "tools")
        ),
        tabItem(tabName = "t_MP",
          HTML('<h4><span style="font-family:Trebuchet MS,Helvetica,sans-serif"><strong><span style="color:#4e5f70">MultiPlots</span></strong></span></h4>'),
          MultiPlotsUI(id = "MP")
        )
      )
    )
  

#### UI ----
ui <- 
  dashboardPage(skin = "blue", 
    header,
    siderbar,
    body
  )
