% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/launch_scX.R
\name{launch_scX}
\alias{launch_scX}
\title{A Shiny app for visualization of single-cell data}
\usage{
launch_scX(
  cseo,
  dataset_name = "scX",
  point.size = 20,
  launch.browser = T,
  host = getOption("shiny.host", "127.0.0.1"),
  port = getOption("shiny.port")
)
}
\arguments{
\item{cseo}{List with a SingleCellExperiment object returned by \code{createSCEobject()}.}

\item{dataset_name}{Optional name for app title.}

\item{point.size}{Point size for plots in app. Defaults to 20.}

\item{launch.browser}{Launch Shiny app in browser or not. Defaults to \code{TRUE}.}
}
\description{
\code{launch_scX()} launches the scX Shiny app for visualizing single-cell data
}
