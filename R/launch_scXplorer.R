#' A Shiny app for visualization of single-cell data
#'
#' `launch_scXplorer()` launches the scXplorer Shiny app for visualizing single-cell data
#'
#' @param cseo List with a SingleCellExperiment object returned by `createSCEobject()`.
#' @param dataset_name Optional name for app title.
#' @param point.size Point size for plots in app. Defaults to 20.
#' @param launch.browser Launch Shiny app in browser or not. Defaults to `TRUE`.
#' @export
launch_scXplorer <- function(cseo, dataset_name='scXplorer', point.size=20, launch.browser = T) {

	options(bitmapType='cairo')
	options(spinner.color="#006272")
  file_path <- system.file("app.R", package = "scXplorer")
  if (!nzchar(file_path)) stop("Shiny app not found")
  ui <- server <- NULL # avoid NOTE about undefined globals
  source(file_path, local = TRUE)
  
  cseo$SCE = cseo$SCE[,cseo$CELLS2KEEP]

  server_env <- environment(server)
  server_env$cseo <- cseo
  server_env$point.size <- point.size
  
  ui_env <- environment(ui)
  ui_env$dataset_name <- dataset_name

  theme_set(theme_bw())
  shiny::shinyApp(ui, server, options = list(launch.browser = launch.browser))
}
# load(file = 'data/Dataset1_paper.Rdata')
#launch_scXplorer(cseo = cseo,
#                 point.size = 20,
#                 dataset_name = 'Dataset 1')
