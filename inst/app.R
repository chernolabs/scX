ui <- function(){
	path_lab <- "https://www.leloir.org.ar/biologia-de-sistemas-integrativa?area=bioinformatica-y-biologia-computacional"
	path_git <- "https://github.com/chernolabs/scX/"
	path_site <- "https://chernolabs.github.io/scX/"
	
	# SideBar ----
	siderbar <- dashboardSidebar(
		width = 250,
		useShinyjs(),
		sidebarMenu(
			id = 'sidebar',
			style = "position: relative; overflow: visible;",
			HTML(
			  paste0(
				"<br>","<a href='",
				path_site,
				"' target='_blank'><img style =
						  'display: block; margin-left: auto; margin-right: auto;'
						  src='www/scX.png' width = '186'></a>",
				# "<br>",
				# paste0(
				  # "<p style = 'text-align: center;'><small><a href='",
				  # path_git,
				  # "' target='_blank'>scX</a></small></p>"
				# ),
				"<br>"
			  )
			),
			tags$style("@import url(https://use.fontawesome.com/releases/v6.2.0/css/all.css);"),
			# menuItem("Global Options", tabName = "globalTab", icon = icon("fa-regular fa-globe")),
			# div( id = 'sidebar_global',
			#      conditionalPanel("input.sidebar === 'globalTab'",
			#                       sliderInput(inputId = "point.size",
			#                                   label = "Point size ScatterPlots",
			#                                   min = 5,max = 50,value = 20),
			#                       fluidRow(column=12, align = "right",
			#                                style='padding-left:12px; padding-right:12px;',
			#                                actionBttn(inputId = "button.config",
			#                                           label = "Apply", 
			#                                           style = "stretch",
			#                                           color = "primary")
			#                       )
			#      )
			# ),
			menuItem("Summary", tabName = "smryTab", icon = icon('bookmark', class = 'fa-solid'),selected = T),
			menuItem("Markers", tabName = "markersTab",
					 icon = icon('location-dot'), startExpanded = F,
					 menuSubItem('Cluster markers', tabName = "markers_cluster",
								 icon = icon('map-location-dot')
					 ),
					 menuSubItem('Find new markers', tabName = "markers_new",
								 icon = icon('magnifying-glass-location')
					 )
			),
			menuItem("Gene Expression", tabName = "g_expTab",
					 icon = icon('chart-simple'), startExpanded = F,
					 menuSubItem('Categories', tabName = "g_exp",
					             icon = icon('signal')
					 ),
					 menuSubItem('Fields', tabName = "g_numeric_exp",
								 icon = icon('audible')
					 ),
					 menuSubItem('Co-expression', tabName = "g_coexp",
								 icon = icon('comments', class = 'fa-regular'))
			),
			menuItem("Differential Expression", tabName = "volcanoTab",
					 icon = icon("chart-column")
			),
			menuItem("Exploratory Data Analysis", tabName = "EDATab",
			         icon = icon("circle-nodes"), startExpanded = F,
			         menuSubItem('Categories', tabName = "clustersTab",
			                     icon = icon('bars-progress')
			         ),
			         menuSubItem('Fields', tabName = "fieldsTab",
			                     icon = icon('clone')
			         )
			),
			menuItem("Visual Tools", tabName = "toolsTab",
					 icon = icon("toolbox"), startExpanded = F,
					 menuSubItem('Violin by Category', tabName = "t_VGL",
								 icon = icon('wrench')
					 ),
					 menuSubItem('MultiPlots', tabName = "t_MP",
								 icon = icon('hammer')
					 )
			)
		)
	)
	
	# Body ----
	body <- dashboardBody(
      useSweetAlert(),
      
      ### Styling ----
      tags$head(
        tags$link(rel="shortcut icon", href="www/scX.ico"),
        tags$script(glue::glue('document.title = {dataset_name}')),
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
		tags$style(HTML('.sidebar .sidebar-menu .treeview-menu>li {margin-left: 10px}')),
        tags$style(
          HTML("hr {border-top: 1px solid #000000;}")
        ),
		# maximum box selection height + scroll (used in multiple gene selection)
		tags$style(".selectize-input {max-height: 100px; overflow-y:auto;}"),
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
	  # tags$footer(align = "left",
			# style = "position:absolute; bottom:0; width:100%; height:20px; z-index:1000;
				# text-align:center; font-size:10px; color:white; padding:1px; background-color:gray;",
			# HTML('<footer> Icon by @irimiodo </footer>')
		# ),
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
                HTML('<h4><span style="font-family:Trebuchet MS,Helvetica,sans-serif"><strong><span style="color:#4e5f70">Expression by Categories</span></strong></span></h4>'),
                ExpressionUI(id="Exp")
        ),
        tabItem(tabName = "g_numeric_exp",
                HTML('<h4><span style="font-family:Trebuchet MS,Helvetica,sans-serif"><strong><span style="color:#4e5f70">Expression by Fields</span></strong></span></h4>'),
                Numeric_ExpressionUI(id="numeric_Exp")
        ),
        tabItem(tabName = "g_coexp",
                HTML('<h4><span style="font-family:Trebuchet MS,Helvetica,sans-serif"><strong><span style="color:#4e5f70">Co-expression</span></strong></span></h4>'),
                COExpUI(id="Co-exp")
        ),
        tabItem(tabName = "volcanoTab",
                HTML('<h4><span style="font-family:Trebuchet MS,Helvetica,sans-serif"><strong><span style="color:#4e5f70">Differential Expression</span></strong></span></h4>'),
                VolcanoUI(id="volcano")
        ),
        tabItem(tabName = "clustersTab",
                HTML('<h4><span style="font-family:Trebuchet MS,Helvetica,sans-serif"><strong><span style="color:#4e5f70">Category Analysis</span></strong></span></h4>'),
                Clusters_UI("cluster")
        ),
        tabItem(tabName = "fieldsTab",
                HTML('<h4><span style="font-family:Trebuchet MS,Helvetica,sans-serif"><strong><span style="color:#4e5f70">Fields Analysis</span></strong></span></h4>'),
                Fields_UI("fields")
        ),
        tabItem(tabName = "t_VGL",
                HTML('<h4><span style="font-family:Trebuchet MS,Helvetica,sans-serif"><strong><span style="color:#4e5f70">Violin by Category</span></strong></span></h4>'),
                VT_UI(id = "tools")
        ),
        tabItem(tabName = "t_MP",
                HTML('<h4><span style="font-family:Trebuchet MS,Helvetica,sans-serif"><strong><span style="color:#4e5f70">MultiPlots</span></strong></span></h4>'),
                MultiPlotsUI(id = "MP")
        )
      )
    )

	# Dashboard UI
	dashboardPage(skin = "blue", 
		header = dashboardHeader(title = HTML(dataset_name), disable = FALSE, titleWidth  = 250),
		sidebar = siderbar,
		body = body)
}

server <- function(input, output, session) {
	QC_Server(id = "qc", sce = cseo$SCE, descriptionText = cseo$text, Npartition = length(cseo$sce.markers))
	markersServer(id="markers", sce=cseo$SCE, sce.markers = cseo$sce.markers, point.size=point.size)
	N_markersServer(id="n_markers", sce=cseo$SCE, point.size=point.size)
	ExpressionServer(id="Exp", sce=cseo$SCE, point.size=point.size)
	Numeric_ExpressionServer(id="numeric_Exp",sce=cseo$SCE,point.size = point.size)
	COExpServer(id="Co-exp", sce=cseo$SCE, point.size = point.size)
	VolcanoServer(id="volcano", sce=cseo$SCE, sce.degs=cseo$sce.degs)
	VT_Server(id = "tools", sce = cseo$SCE)
	MultiPlotsServer(id = "MP", sce=cseo$SCE)
	Clusters_Server(id ="cluster", sce=cseo$SCE)
	Fields_Server("fields",sce = cseo$SCE)                                                                                                                                                                                                     
}
