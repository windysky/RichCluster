#
# Shiny web application for integrative enrichment analysis and visualization
#

library(shiny)
library(Rcpp)
library(shinydashboard)
library(DT)
library(plotly)


ui <- function(request) {
  dashboardPage(
    dashboardHeader(title = "RichCluster"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Home", icon=icon("house"), tabName="home_tab"),
        menuItem("Cluster", icon=icon("th"), tabName="cluster_tab"),
        menuItem("Visualize", icon=icon("eye"), tabName="visualize_tab")
      )
    ),
    dashboardBody(
    tabItems(
      homeTabUI("home", tabName="home_tab"),
      clusterTabUI("cluster", tabName="cluster_tab"),
      visualizeTabUI("visualize", tabName="visualize_tab")
    ),
    tags$head(
      tags$style(
        HTML('
             .content-wrapper { overflow: auto; }
             .dataTables_wrapper { overflow-x: scroll; }
            ')
        )
      )
    )
  )
}

server <- function(input, output, session) {

  richsets <- reactiveValues()
  richnames <- reactiveValues(labels=NULL)

  homeTabServer("home")
  clusterTabServer("cluster", richsets, richnames)
  visualizeTabServer("visualize")

}
# Run the application
shinyApp(ui=ui, server=server)



